#include "Processor.h"
#include <stdexcept>

#include <cassert>


using namespace std;
using namespace ramulator;

Processor::Processor(const Config& configs,
    vector<string> trace_list,
    function<bool(Request)> send_memory,
    MemoryBase& memory)
    : early_exit(configs.is_early_exit()),
    no_core_caches(!configs.has_core_caches()),
    no_shared_cache(!configs.has_l3_cache()),
    cachesys(new CacheSystem(configs, send_memory)),
    llc(l3_size, l3_assoc, l3_blocksz,
         mshr_per_bank * trace_list.size(),
         Cache::Level::L3, cachesys), memory(memory){

  //set initial parameters
  assert(cachesys != nullptr);
  pim_mode_enabled = configs.pim_mode_enabled();
  trace.pim_mode_enabled = pim_mode_enabled;
  pim_cores = configs.get_num_pim_cores();

  if(configs.get_trace_format() == Config::Format::PISA)
      pisa_trace = true;
  else if (configs.get_trace_format() == Config::Format::ZSIM)
      zsim_trace = true;
  cycle_time = configs.get_cpu_tick()/1000.0;

  //create cores
  if (no_shared_cache) {
      if(configs.pim_mode_enabled()){
          cout << "Number of Cores: " << pim_cores << endl;
          for (int i = 0 ; i < pim_cores ; ++i) {
            cores.emplace_back(new Core(configs, i, send_memory, nullptr, cachesys, memory));
          }
      }
  }
  else {
    for (int i = 0 ; i < pim_cores ; i++) {
      cores.emplace_back(new Core(configs, i,
          std::bind(&Cache::send, &llc, std::placeholders::_1),
          &llc, cachesys, memory));
    }
  }

  //check if single or multiple traces
  if(!configs.get_split_trace()){
    cout << "Single trace \n";
    trace.init_trace(trace_list[0]);
    int tracenum = trace_list.size();
    assert(tracenum > 0);
    printf("tracenum: %d\n", tracenum);
    for (int i = 0 ; i < tracenum ; ++i) {
      printf("trace_list[%d]: %s\n", i, trace_list[i].c_str());
    }
  }
  else{
      cout << "Split Trace \n";
      for (int i = 0 ; i < pim_cores ; ++i) {
        cores[i]->load_trace(trace_list[0]);
        cores[i]->split_trace = true;
        cores[i] -> get_first_request();
      }
  }

  //bind cores to memory
  if(configs.pim_mode_enabled() == 1){
      ipcs.resize(pim_cores);
      for (int i = 0 ; i < pim_cores ; ++i) {
        cores[i]->callback = std::bind(&Processor::receive, this,
            placeholders::_1);
        ipcs[i] = -1;
      }
  }

  //if pim_mode and pin trace and single trace
  if(pim_mode_enabled && !pisa_trace && !zsim_trace && !configs.get_split_trace()){
    long bubble_cnt;
    long req_addr = -1;
    Request::Type req_type;
    bool more_reqs = trace.get_unfiltered_request(bubble_cnt, req_addr, req_type);

    while(more_reqs){
        MemoryBase *ptr_mem = &memory;
        Memory<HMC,Controller>* ptr = dynamic_cast<Memory<HMC,Controller>*>(ptr_mem);
        if(ptr != NULL){
            long tmp = req_addr;
            ptr -> clear_higher_bits(tmp, ptr->max_address-1ll);
            ptr -> clear_lower_bits(tmp, ptr -> tx_bits);
            int max_block_col_bits =  ptr->spec->maxblock_entry.flit_num_bits - ptr->tx_bits;
            ptr->slice_lower_bits(tmp, max_block_col_bits);
            int vault_target = ptr->slice_lower_bits(tmp, ptr->addr_bits[int(HMC::Level::Vault)]);

            Core* core = cores[vault_target % pim_cores].get();
            core -> rqst_queue.bubble_cnt.push_back(bubble_cnt);
            core -> rqst_queue.req_addr.push_back(req_addr);
            core -> rqst_queue.req_type.push_back(req_type);
            core -> rqst_queue.numberInstructionsInQueue++;
        }
        else{
            cout << "Bad conversion \n";
        }
        more_reqs = trace.get_unfiltered_request(bubble_cnt, req_addr, req_type);
    }
    for(int i = 0; i < cores.size(); i++){
          Core* core = cores[i].get();
          cout << "Core " << i << " got " << core->rqst_queue.numberInstructionsInQueue << endl;
          if(!core -> rqst_queue.is_empty()){
              core -> bubble_cnt = core -> rqst_queue.bubble_cnt.front();
              core -> req_addr = core -> rqst_queue.req_addr.front();
              core -> req_type = core -> rqst_queue.req_type.front();
              core -> more_reqs = true;
              core->rqst_queue.pop_front();

          }
          else{
              core -> more_reqs = false;
          }
    }
    cout << endl;

  }
  //if pim mode and pisa trace and single trace
  else if(pim_mode_enabled && pisa_trace && !zsim_trace && !configs.get_split_trace()){
    long bubble_cnt;
    long req_addr = -1;
    unsigned size = 0;
    Request::Type req_type;
    long totalReads = 0;
    long totalWrites = 0;
    bool more_reqs = trace.get_pisa_request(bubble_cnt, req_addr, req_type, size);
    (req_type == Request::Type::READ) ? totalReads++ : totalWrites++;

    vector<int> round_robin(32, 0);
    int number_cores_per_vault = pim_cores/32;
    while(more_reqs){
        MemoryBase *ptr_mem = &memory;
        Memory<HMC,Controller>* ptr = dynamic_cast<Memory<HMC,Controller>*>(ptr_mem);
        if(ptr != NULL){
            long tmp = req_addr;
            ptr -> clear_higher_bits(tmp, ptr->max_address-1ll);
            ptr -> clear_lower_bits(tmp, ptr -> tx_bits);
            int max_block_col_bits =  ptr->spec->maxblock_entry.flit_num_bits - ptr->tx_bits;
            ptr->slice_lower_bits(tmp, max_block_col_bits);
            int vault_target = ptr->slice_lower_bits(tmp, ptr->addr_bits[int(HMC::Level::Vault)]);

            if(pim_cores <= 32){
                Core* core = cores[vault_target % pim_cores].get();
                core -> rqst_queue.bubble_cnt.push_back(bubble_cnt);
                core -> rqst_queue.req_addr.push_back(req_addr);
                core -> rqst_queue.req_type.push_back(req_type);
                core -> rqst_queue.numberInstructionsInQueue++;
            }
            else{
                int core_to_schedule = 32*round_robin[vault_target] + vault_target;
                round_robin[vault_target] = (round_robin[vault_target] + 1) % number_cores_per_vault;
                if(core_to_schedule < pim_cores){
                    Core* core = cores[core_to_schedule].get();
                    core -> rqst_queue.bubble_cnt.push_back(bubble_cnt);
                    core -> rqst_queue.req_addr.push_back(req_addr);
                    core -> rqst_queue.req_type.push_back(req_type);
                    core -> rqst_queue.numberInstructionsInQueue++;
                }
            }
        }
        else{
            cout << "Bad conversion \n";
        }
        more_reqs = trace.get_pisa_request(bubble_cnt, req_addr, req_type,size);
        (req_type == Request::Type::READ) ? totalReads++ : totalWrites++;
    }
    cout << "Total Reads: " << totalReads << " Total Writes: " << totalWrites << endl;
    for(int i = 0; i < cores.size(); i++){
          Core* core = cores[i].get();
          cout << "Core " << i << " got " << core->rqst_queue.numberInstructionsInQueue << endl;
          if(!core -> rqst_queue.is_empty()){
              core -> bubble_cnt = core -> rqst_queue.bubble_cnt.front();
              core -> req_addr = core -> rqst_queue.req_addr.front();
              core -> req_type = core -> rqst_queue.req_type.front();
              core -> more_reqs = true;
              core->rqst_queue.pop_front();

          }
          else{
              core -> more_reqs = false;
          }
    }
    cout << endl;
  }
  //if pim mode and zsim trace and single trace
  else if(pim_mode_enabled && !pisa_trace && zsim_trace && !configs.get_split_trace()){
    long bubble_cnt;
    long req_addr = -1;
    unsigned cpu_id = 0;
    Request::Type req_type;
    long totalReads = 0;
    long totalWrites = 0;
    if(configs.get_disable_per_scheduling() == false){
        cout << "Perfect scheduling enabled \n";
    }
    else{
        cout << "Perfect scheduling disabled \n";
    }
    bool more_reqs = trace.get_zsim_request(bubble_cnt, req_addr, req_type, cpu_id);
    (req_type == Request::Type::READ) ? totalReads++ : totalWrites++;
    vector<int> round_robin(32, 0);
    int number_cores_per_vault = pim_cores/32;
    vector<int> offset_counter(256,0);
    while(more_reqs){
        if(req_addr != -1){
            if(configs.get_disable_per_scheduling() == false){
                MemoryBase *ptr_mem = &memory;
                Memory<HMC,Controller>* ptr = dynamic_cast<Memory<HMC,Controller>*>(ptr_mem);
                if(ptr != NULL){
                    req_addr = ptr->page_allocator(req_addr, 0);
                    long tmp = req_addr;
                    ptr -> clear_higher_bits(tmp, ptr->max_address-1ll);
                    ptr -> clear_lower_bits(tmp, ptr -> tx_bits);
                    int max_block_col_bits =  ptr->spec->maxblock_entry.flit_num_bits - ptr->tx_bits;
                    ptr->slice_lower_bits(tmp, max_block_col_bits);
                    int vault_target = ptr->slice_lower_bits(tmp, ptr->addr_bits[int(HMC::Level::Vault)]);

                    if(pim_cores <= 32){
                        Core* core = cores[vault_target % pim_cores].get();
                        core -> rqst_queue.bubble_cnt.push_back(bubble_cnt);
                        core -> rqst_queue.req_addr.push_back(req_addr);
                        core -> rqst_queue.req_type.push_back(req_type);
                        core -> rqst_queue.numberInstructionsInQueue++;
                    }
                    else{
                        int core_to_schedule = 32*round_robin[vault_target] + vault_target;
                        round_robin[vault_target] = (round_robin[vault_target] + 1) % number_cores_per_vault;
                        if(core_to_schedule < pim_cores){
                            Core* core = cores[core_to_schedule].get();
                            core -> rqst_queue.bubble_cnt.push_back(bubble_cnt);
                            core -> rqst_queue.req_addr.push_back(req_addr);
                            core -> rqst_queue.req_type.push_back(req_type);
                            core -> rqst_queue.numberInstructionsInQueue++;
                        }
                    }
                }
                else{
                    cout << "Bad conversion \n";
                }
            }
            else{
                if(pim_cores <= 256){
                    req_addr = memory.page_allocator(req_addr, cpu_id);
                    Core* core = cores[cpu_id % pim_cores].get();
                    core -> rqst_queue.bubble_cnt.push_back(bubble_cnt);
                    core -> rqst_queue.req_addr.push_back(req_addr);
                    core -> rqst_queue.req_type.push_back(req_type);
                    core -> rqst_queue.numberInstructionsInQueue++;
                }
                else{
                    int new_cpu_id = cpu_id + offset_counter[cpu_id]*256;
                    offset_counter[cpu_id]++;
                    req_addr = memory.page_allocator(req_addr, new_cpu_id);
                    Core* core = cores[new_cpu_id % pim_cores].get();
                    core -> rqst_queue.bubble_cnt.push_back(bubble_cnt);
                    core -> rqst_queue.req_addr.push_back(req_addr);
                    core -> rqst_queue.req_type.push_back(req_type);
                    core -> rqst_queue.numberInstructionsInQueue++;
                }
            }
        }
        more_reqs = trace.get_zsim_request(bubble_cnt, req_addr, req_type,cpu_id);
        (req_type == Request::Type::READ) ? totalReads++ : totalWrites++;
    }

    cout << "Total Reads: " << totalReads << " Total Writes: " << totalWrites << endl;
    for(int i = 0; i < cores.size(); i++){
          Core* core = cores[i].get();
          cout << "Core " << i << " got " << core->rqst_queue.numberInstructionsInQueue << endl;
          if(!core -> rqst_queue.is_empty()){
              core -> bubble_cnt = core -> rqst_queue.bubble_cnt.front();
              core -> req_addr = core -> rqst_queue.req_addr.front();
              core -> req_type = core -> rqst_queue.req_type.front();
              core -> more_reqs = true;
              core->rqst_queue.pop_front();

          }
          else{
              core -> more_reqs = false;
          }
    }
    cout << endl;
  }

  expected_limit_insts = configs.get_expected_limit_insts();


  //done with configuration, set stats
  // regStats
  cpu_cycles.name("cpu_cycles")
            .desc("cpu cycle number")
            .precision(0)
            ;

  // regStats
  general_ipc.name("ipc")
            .desc("final ipc number")
            .precision(6)
            ;
// regStats
  total_time.name("total_time")
            .desc("Total Time (ns)")
            .precision(6)
            ;
  // regStats
  total_cpu_instructions.name("cpu_instructions")
            .desc("total cpu instructions number")
            .precision(0)
            ;

  // regStats
  total_cache_misses.name("total_cache_misses_pim")
            .desc("Total cache misses at pim")
            .precision(0)
            ;

  // regStats
  total_cache_hits.name("total_cache_hits_pim")
            .desc("Total cache hits at pim")
            .precision(0)
            ;
  // regStats
  total_idle_cycles.name("Total idle cycles")
            .desc("Total idle cycles due to full window.")
            .precision(0)
            ;
  // regStats
  average_idle_cycles.name("Average idle cycles")
            .desc("Average idle cycles due to full window")
            .precision(0)
            ;

  general_ipc = 0.0;
  total_cpu_instructions = 0;
  cpu_cycles = 0;
  average_idle_cycles = 0;
  total_idle_cycles = 0;
  total_cache_misses = 0;
  total_cache_hits = 0;
}

void Processor::tick() {
  cpu_cycles++;

   if((int(cpu_cycles.value()) % 1000000) == 0)
      printf("CPU heartbeat, cycles: %d \n", (int(cpu_cycles.value())));
  /*    for (unsigned int i = 0 ; i < cores.size() ; ++i) {
        Core* core = cores[i].get();
        if(core->finished())
            cout << "Core " << i << " has finished \n";
        else{
            //if(core->rqst_queue.numberInstructionsInQueue == core->prevnumberInstructionsInQueue){
            //    exit(1);
            //}
            
            core->prevnumberInstructionsInQueue = core->rqst_queue.numberInstructionsInQueue;
            cout << "Core " << i << " hasn't finished  - "
                                    "number of instructions: " << core->rqst_queue.numberInstructionsInQueue << " \n";

        }
      }
      cout << endl;
  }*/

  if (!(no_core_caches && no_shared_cache)) {
    //cout << "Ticking Cache \n";
    cachesys->tick();
  }

  for (unsigned int i = 0 ; i < cores.size() ; i++) {
    Core* core = cores[i].get();
    core->tick();
  }
}

void Processor::receive(Request& req) {
  if (!no_shared_cache) {
    llc.callback(req);
  } else if (!cores[0]->no_core_caches) {
    // Assume all cores have caches or don't have caches
    // at the same time.
    for (unsigned int i = 0 ; i < cores.size() ; ++i) {
      Core* core = cores[i].get();
      core->caches[0]->callback(req);
    }
  }

    //cout << "Core " << req.coreid << " receive response " << req.initial_addr << "\n";

    Core* core = cores[req.coreid].get();
    core->receive(req);

   //for (unsigned int i = 0 ; i < cores.size() ; ++i) {
   //  Core* core = cores[i].get();
   // core->receive(req);
  //    }
}

void Processor::calc_stats(){
for (unsigned int i = 0 ; i < cores.size(); ++i){
      if (ipcs[i] < 0) {
        ipcs[i] = cores[i]->calc_ipc();
        ipc += ipcs[i];
        total_retired  += cores[i]->retired;
        total_instructions += cores[i]-> cpu_inst.value();
        total_idle_cycles += cores[i] -> idle_cycles.value();
        total_cache_hits += cores[i] -> cache_hits;
        total_cache_misses += cores[i] -> cache_misses;
      }
      else{
          return;
      }
    }

    ipc = total_instructions/cpu_cycles.value();
    average_idle_cycles = total_idle_cycles.value()/cores.size();

    cout << endl;
    cout << "-> retired: " << total_retired << endl;
    cout << "-> cycles: " << cpu_cycles.value() << endl;
    cout << "-> ipc: " << ipc << endl;
    cout << "-> total instructions: " << total_instructions << endl;

    general_ipc = ipc;
    total_cpu_instructions = total_instructions;
    if(pim_mode_enabled)
        total_time = total_instructions*(1/ipc)*0.8;
    else
        total_time = total_instructions*(1/ipc)*cycle_time;
    cout << "-> total time: " << total_time.value() << "ns" <<endl;

}
bool Processor::finished() {
  if (early_exit) {
    for (unsigned int i = 0 ; i < cores.size(); ++i) {
      if (cores[i]->finished()) {
        for (unsigned int j = 0 ; j < cores.size() ; ++j) {
          ipc += cores[j]->calc_ipc();
        }
        return true;
      }
    }
    return false;
  }
  else {
    for (unsigned int i = 0 ; i < cores.size(); ++i) {
      if (!cores[i]->finished()) {
        return false;
      }
    }

    calc_stats();
    return true;
  }
}

bool Processor::has_reached_limit() {

   long total_instr = 0;
   for (unsigned int i = 0 ; i < cores.size() ; ++i) {
      total_instr += long(cores[i]->cpu_inst.value());
    }

    if(total_instr >= expected_limit_insts) return true;



  for (unsigned int i = 0 ; i < cores.size() ; ++i) {
    if (!cores[i]->has_reached_limit()) {
      return false;
    }
  }
  calc_stats();
  return true;
}

Core::Core(const Config& configs, int coreid, function<bool(Request)> send_next,
    Cache* llc, std::shared_ptr<CacheSystem> cachesys, MemoryBase& memory)
    : id(coreid), no_core_caches(!configs.has_core_caches()),
    no_shared_cache(!configs.has_l3_cache()),
    llc(llc), memory(memory){


  cpu_type = configs.get_cpu_type();
  inFlightMemoryAccess = 0;
  //window.depth = configs.get_window_depth();
  /*if(configs.get_pim_cacheline_size() < ways){
      ways = configs.get_pim_cacheline_size();
      number_cache_lines = configs.get_pim_cacheline_size();
  }
  else{
      number_cache_lines = configs.get_pim_cacheline_size()/ways;
  }

  pim_private_cache.resize(ways);
  for(int i=0; i < ways; i++){
    pim_private_cache[i].valid.resize(number_cache_lines, false);
    pim_private_cache[i].lru_counter.resize(number_cache_lines, i);
    pim_private_cache[i].tag.resize(number_cache_lines);
 }*/

  if (no_core_caches) {
    send = send_next;
  } else {
    // L2 caches[0]
    caches.emplace_back(new Cache(
        l2_size, l2_assoc, l2_blocksz, l2_mshr_num,
        Cache::Level::L2, cachesys));
    // L1 caches[1]
    caches.emplace_back(new Cache(
        l1_size, l1_assoc, l1_blocksz, l1_mshr_num,
        Cache::Level::L1, cachesys));
    send = bind(&Cache::send, caches[1].get(), placeholders::_1);
    if (llc != nullptr) {
      caches[0]->concatlower(llc);
    }
    caches[1]->concatlower(caches[0].get());
  }

  
  // set expected limit instruction for calculating weighted speedup
  expected_limit_insts = configs.get_expected_limit_insts();

  //cout << "expected limit insts: " << expected_limit_insts << endl;
  // regStats
  record_cycs.name("record_cycs_core_" + to_string(id))
             .desc("Record cycle number for calculating weighted speedup. (Only valid when expected limit instruction number is non zero in config file.)")
             .precision(0)
             ;

  record_insts.name("record_insts_core_" + to_string(id))
              .desc("Retired instruction number when record cycle number. (Only valid when expected limit instruction number is non zero in config file.)")
              .precision(0)
              ;



  memory_access_cycles.name("memory_access_cycles_core_" + to_string(id))
                      .desc("memory access cycles in memory time domain")
                      .precision(0)
                      ;
  memory_access_cycles = 0;

  cpu_inst.name("cpu_instructions_core_" + to_string(id))
          .desc("cpu instruction number")
          .precision(0)
          ;
  cpu_inst = 0;

  idle_cycles.name("idle_cycles_core_" + to_string(id))
          .desc("idle cycles due to windows full")
          .precision(0)
          ;
  idle_cycles = 0;

    memory_inst.name("memory_instructions_core_"+ to_string(id))
                    .desc("memory instruction number")
                    .precision(0)
          ;

    non_memory_inst.name("non_memory_instructions_core_"+ to_string(id))
                    .desc("non memory instructions")
                    .precision(0);


     cpu_inst = 0;
     memory_inst = 0;
     non_memory_inst = 0;

}

void Core::load_trace(string trace_base_name){
    cout << "Core " << id << " trying to load trace " << trace_base_name+"."+std::to_string(id) << endl;
    trace_per_core.init_trace(trace_base_name+"."+std::to_string(id));
}

void Core::get_first_request(){
    unsigned int cpu_id;
    more_reqs = trace_per_core.get_zsim_request(bubble_cnt, req_addr, req_type,cpu_id);
    req_addr = memory.page_allocator(req_addr, id);
}
double Core::calc_ipc()
{
    printf("[%d]retired: %ld, clk, %ld\n", id, retired, clk);
    return (double) retired / clk;
}


void Core::tick_outOrder(){
    retired += window.retire();
    //if(retired == 0 ) {cout << "Core " << id << " retired 0 instructions \n";}
    if (expected_limit_insts == 0 && !more_reqs){
        return;
    }

    if(lock_core){return;}

    // bubbles (non-memory operations)
    int inserted = 0;
    while (bubble_cnt > 0) {
        if (inserted == window.ipc){idle_cycles++; return;}
        if (window.is_full()){ return;}

        window.insert(true, -1);
        inserted++;
        bubble_cnt--;
        cpu_inst++;
        non_memory_inst++;
        if (long(cpu_inst.value()) == expected_limit_insts && !reached_limit) {
          record_cycs = clk;
          record_insts = long(cpu_inst.value());
          memory.record_core(id);
          reached_limit = true;
        }
    }

    if (req_type == Request::Type::READ) {
        // read request
        if (inserted == window.ipc) {idle_cycles++; return;}
        if (window.is_full()){ return;}

       Request req(req_addr, req_type, callback, id);
       if (!send(req)) return;
      // cout << "Core " << id << " sent " << req_addr << endl;
       window.insert(false, req_addr);
       cpu_inst++;
       memory_inst++;
    }
    else if(req_type == Request::Type::WRITE){
        // write request
        assert(req_type == Request::Type::WRITE);

        Request req(req_addr, req_type, callback, id);
        if (!send(req)) return;
        cpu_inst++;
        memory_inst++;

    }
    else if(req_type == Request::Type::INSTRUCTION){
        if (inserted == window.ipc) {idle_cycles++; return;}

        Request req(req_addr, Request::Type::READ, callback, id);
        req.instruction_request = true;
        if (!send(req)) return;
        //lock_core = true;
    }

    //get next request
    if(split_trace){
        unsigned int cpu_id;
        more_reqs = trace_per_core.get_zsim_request(bubble_cnt, req_addr, req_type,cpu_id);
        req_addr = memory.page_allocator(req_addr, id);
    }
    else{
        if(!rqst_queue.is_empty()){
            more_reqs = true;
            bubble_cnt = rqst_queue.bubble_cnt.front();
           // cout << "bubble count: " << bubble_cnt << endl;
            req_addr = rqst_queue.req_addr.front();
            req_type = rqst_queue.req_type.front();
            rqst_queue.pop_front();
        }
        else{
            more_reqs = false;
        }
    }

    if (!more_reqs) {
       if (!reached_limit) { // if the length of this trace is shorter than expected length, then record it when the whole trace finishes, and set reached_limit to true.
         record_cycs = clk;
         record_insts = long(cpu_inst.value());
         memory.record_core(id);
         reached_limit = true;
       }
     }
}

void Core::tick_inOrder(){

    if (expected_limit_insts == 0 && !more_reqs){
        return;
    }

    if(inFlightMemoryAccess >= 1){
        return;
    }

    if(lock_core){cout << "Core locked \n"; return; }


    // bubbles (non-memory operations)
    int inserted = 0;
    while (bubble_cnt > 0) {
        if (inserted == window.ipc){return;}
        inserted++;
        bubble_cnt--;
        cpu_inst++;
        retired++;

        if (long(cpu_inst.value()) == expected_limit_insts && !reached_limit) {
          record_cycs = clk;
          record_insts = long(cpu_inst.value());
          memory.record_core(id);
          reached_limit = true;
        }
    }


    if(check_cache_hit(req_addr) && (req_type == Request::Type::READ)){
        if (inserted == window.ipc) {return;}
        cpu_inst++;
    }
    else if (req_type == Request::Type::READ) {
        // read request
       if (inserted == window.ipc) {return;}

       Request req(req_addr, req_type, callback, id);

       if (!send(req)) {return;}
       inFlightMemoryAccess++;
       cpu_inst++;
       if (long(cpu_inst.value()) == expected_limit_insts && !reached_limit) {
         record_cycs = clk;
         record_insts = long(cpu_inst.value());
         memory.record_core(id);
         reached_limit = true;
       }
    }
    else if(req_type == Request::Type::WRITE){
        // write request
        assert(req_type == Request::Type::WRITE);
        Request req(req_addr, req_type, callback, id);
        if (!send(req)){return;}
        cpu_inst++;
    }
    else if(req_type == Request::Type::INSTRUCTION){
        cout << "Here \n";
        if (inserted == window.ipc) return;

        Request req(req_addr, Request::Type::READ, callback, id);
        req.instruction_request = true;

        if (!send(req)) return;
        inFlightMemoryAccess++;
        lock_core = true;
    }

    if(!rqst_queue.is_empty()){
        more_reqs = true;
        bubble_cnt = rqst_queue.bubble_cnt.front();
        req_addr = rqst_queue.req_addr.front();
        req_type = rqst_queue.req_type.front();
        rqst_queue.pop_front();
    }
    else{
        more_reqs = false;
    }


    if (!more_reqs) {
       if (!reached_limit) { // if the length of this trace is shorter than expected length, then record it when the whole trace finishes, and set reached_limit to true.
         record_cycs = clk;
         record_insts = long(cpu_inst.value());
         memory.record_core(id);
         reached_limit = true;
       }
     }


}

void Core::tick()
{
    clk++;
    if(cpu_type == "inOrder")
        tick_inOrder();
    else if(cpu_type == "outOrder")
        tick_outOrder();
    else
        cout << "Something is wrong \n";
}

bool Core::finished()
{
    return !more_reqs; 
}

void Core::update_lru(unsigned index_number, unsigned way_number ){
    unsigned lru_value = pim_private_cache[way_number].lru_counter[index_number];

    for(int i = 0; i < ways; i++){
        if(way_number == i){
            pim_private_cache[i].lru_counter[index_number] = ways -1;
        }
        else{
            if(pim_private_cache[i].lru_counter[index_number] > lru_value )
                pim_private_cache[i].lru_counter[index_number]--;
        }
    }
}

bool Core::check_cache_hit(long addr){
    unsigned offset = 8;
    int index_number = addr >> offset;
    int index_size = log2(number_cache_lines);
    long lbits = (1 << index_size) - 1;
    index_number &= lbits;

    if(index_number > number_cache_lines){
        cout << "@Core::check_cache_hit - invalid index number \n";
        exit(1);
    }
    for(int i = 0; i < ways; i++){
        if(((addr >> (offset+index_size)) == (pim_private_cache[i].tag[index_number] >> (offset+index_size))) && (pim_private_cache[i].valid[index_number])){
            update_lru(index_number,i);
            cache_hits++;
            //cout << "cache hit \n";
            return true;
        }
    }

    cache_misses++;
    return false;
}

bool Core::has_reached_limit() {
  return reached_limit;
}

void Core::insert_cache_line(long addr){
    int insert_line;
    unsigned offset = 8;
    int index_number = addr >> offset;
    int index_size = log2(number_cache_lines);
    long lbits = (1 << index_size) - 1;
    index_number &= lbits;

    if(index_number > number_cache_lines){
        cout << "@Core::insert_cache_line - invalid index number \n";
        exit(1);
    }
    for(int i = 0; ways; i++){
        if(pim_private_cache[i].lru_counter[index_number] == 0 ){
            insert_line = i;
            break;
        }
    }

    for(int i = 0; i < ways; i++){
        if(i == insert_line ){
            pim_private_cache[i].tag[index_number] = addr;
            pim_private_cache[i].lru_counter[index_number] = ways -1;
            pim_private_cache[i].valid[index_number] = true;
        }
        else{
            pim_private_cache[i].lru_counter[index_number]--;
        }
    }
}
void Core::receive(Request& req){
  //  if(pim_mode_enabled){
        //if(req.coreid == id){
         //   MemoryBase *ptr_mem = &memory;
         //   Memory<HMC,Controller>* ptr = dynamic_cast<Memory<HMC,Controller>*>(ptr_mem);
         //   if(ptr != NULL){
         //       ptr -> read_latency_sum += req.depart_hmc - req.arrive_hmc;
         //       ptr -> request_packet_latency_sum += req.arrive - req.arrive_hmc;
         //       ptr -> response_packet_latency_sum += req.depart_hmc - req.depart;
         //   }
         //   else{
         //       cerr << "Bad Conversion! \n";
         //   }

           //cout << "Got request back, addr: " << req.initial_addr << " Core " << id << endl;
           //insert_cache_line(req.initial_addr);
           //if(cpu_type == "outOrder" && req.type == Request::Type::READ){
        //   cout << "@ Core receive \n";
                window.set_ready(req.addr, ~(l1_blocksz - 1l));

           //}
           //else if(cpu_type == "inOrder"){
           //     inFlightMemoryAccess--;
           //    // cout << "Freeing inflight memory access from Core : " << id << endl;
          // }

           //if(req.instruction_request)
           //    lock_core = false;

           if (req.arrive != -1 && req.depart > last) {
                memory_access_cycles += (req.depart - max(last, req.arrive));
                last = req.depart;
            }
       // }
    //}
  //  else{
   //     cout <<"here wrong \n";
   //     window.set_ready(req.initial_addr, ~(l1_blocksz - 1l));
   //     if (req.arrive != -1 && req.depart > last) {
   //         memory_access_cycles += (req.depart - max(last, req.arrive));
  //          last = req.depart;
  //      }
   //}
}

bool Window::is_full()
{
    return load == depth;
}

bool Window::is_empty()
{
    return load == 0;
}


void Window::insert(bool ready, long addr)
{
    assert(load <= depth);

    ready_list.at(head) = ready;
    addr_list.at(head) = addr;

    head = (head + 1) % depth;
    load++;
}


long Window::retire()
{
    assert(load <= depth);

    if (load == 0) return 0;

    int retired = 0;
    while (load > 0 && retired < ipc) {
        if (!ready_list.at(tail))
            break;

        tail = (tail + 1) % depth;
        load--;
        retired++;
    }

    return retired;
}


void Window::set_ready(long addr, int mask)
{
    if (load == 0) return;

    for (int i = 0; i < load; i++) {
        int index = (tail + i) % depth;
        if ((addr_list.at(index) & mask) != (addr & mask))
            continue;
        ready_list.at(index) = true;
    }
}

Trace::Trace(const string& trace_fname) : file(trace_fname), trace_name(trace_fname){

#ifndef MAX_NUMBER_CORES
    #define MAX_NUMBER_CORES 16
#endif

    instructions.resize(MAX_NUMBER_CORES);
    for(int i = 0; i < MAX_NUMBER_CORES; i++) instructions[i] = 0;
    if (!file.good()) {
        std::cerr << "Bad trace file: " << trace_fname << std::endl;
        exit(1);
    }
}

bool Trace::init_trace(const string& trace_fname){

#ifndef MAX_NUMBER_CORES
    #define MAX_NUMBER_CORES 16
#endif

    instructions.resize(MAX_NUMBER_CORES);
    for(int i = 0; i < MAX_NUMBER_CORES; i++) instructions[i] = 0;
    trace_name = trace_fname;
    file.open(trace_name);

    if (!file.good()) {
        std::cerr << "Bad trace file: " << trace_fname << std::endl;
        return false;
    }
    cout << "Trace opended: " << trace_fname << endl;
    return true;
}

bool Trace::get_unfiltered_request(long& bubble_cnt, long& req_addr, Request::Type& req_type)
{
    string line;
    getline(file, line);
    if (file.eof()) {
      if(!pim_mode_enabled){
        file.clear();
        file.seekg(0, file.beg);
      }
       return false;
    }

    size_t pos, end;
    bubble_cnt = std::stoul(line, &pos, 10);
    pos = line.find_first_not_of(' ', pos+1);
    req_addr = std::stoul(line.substr(pos), &end, 0);

    pos = line.find_first_not_of(' ', pos+end);

    if (pos == string::npos || line.substr(pos)[0] == 'R')
        req_type = Request::Type::READ;
    else if (line.substr(pos)[0] == 'W')
        req_type = Request::Type::WRITE;
    else{
        cout << "line: " << line << endl;
        cout << line.substr(pos)[0] << endl;
        assert(false);
    }
    return true;
}

bool Trace::get_filtered_request(long& bubble_cnt, long& req_addr, Request::Type& req_type)
{
    static bool has_write = false;
    static long write_addr;
    static int line_num = 0;

    if (has_write){
        bubble_cnt = 0;
        req_addr = write_addr;
        req_type = Request::Type::WRITE;
        has_write = false;
        return true;
    }
    string line;
    getline(file, line);
    line_num ++;
    if (file.eof() || line.size() == 0) {
        if(!pim_mode_enabled){
            file.clear();
            file.seekg(0, file.beg);
        }
        has_write = false;
        line_num = 0;

        if(expected_limit_insts == 0) {
            has_write = false;
            return false;
        }
        else { // starting over the input trace file
            getline(file, line);
            line_num++;
        }
    }

    size_t pos, end;
    bubble_cnt = std::stoul(line, &pos, 10);

    pos = line.find_first_not_of(' ', pos+1);
    req_addr = stoul(line.substr(pos), &end, 0);


    req_type = Request::Type::READ;

    pos = line.find_first_not_of(' ', pos+end);
    if (pos != string::npos){
        has_write = true;
        write_addr = stoul(line.substr(pos), NULL, 0);
    }
    return true;
}

bool Trace::get_dramtrace_request(long& req_addr, Request::Type& req_type)
{
    string line;
    getline(file, line);
    if (file.eof()) {
        return false;
    }
    size_t pos;
    req_addr = std::stoul(line, &pos, 16);

    pos = line.find_first_not_of(' ', pos+1);

    if (pos == string::npos || line.substr(pos)[0] == 'R')
        req_type = Request::Type::READ;
    else if (line.substr(pos)[0] == 'W')
        req_type = Request::Type::WRITE;
    else assert(false);
    return true;
}

std::vector<std::string> Trace::split_str(std::string to_split){
    std::vector<std::string> array;
   // std::size_t pos = 0, found;

    std::istringstream iss(to_split);
        for(std::string to_split; iss >> to_split; )
            array.push_back(to_split);

    /*while((found = to_split.find_first_of(' ', pos)) != std::string::npos) {
        array.push_back(to_split.substr(pos, found - pos));
        pos = found+1;
    }
    array.push_back(to_split.substr(pos));*/
    return array;
}


bool Trace::get_pisa_request(long& bubble_cnt, long& req_addr, Request::Type& req_type, unsigned& size){
    string linebuffer;

    if(getline(file, linebuffer)){
        if (linebuffer.length() == 0) return false;

        vector<string> split = split_str(linebuffer);
        line l;
        l.threadID = std::stoi(split[0],nullptr,10);
        l.processorID = std::stoi(split[1],nullptr,10);
        l.cycleNum = std::stoi(split[2],nullptr,10);

        int nInstructions = l.cycleNum - instructions[l.processorID];
        instructions[l.processorID] = l.cycleNum;
        l.cycleNum = nInstructions;

        l.type = split[3];
        l.addr  = std::stoul(split[4], nullptr, 16);
        l.size = std::stoi(split[5],nullptr,10);

        bubble_cnt = l.cycleNum;
        req_addr = l.addr;
        (l.type == "L")? req_type = Request::Type::READ : req_type = Request::Type::WRITE;
        size = l.size;

        return true;
      }
    else{
        return false;
    }
}

bool Trace::get_zsim_request(long& bubble_cnt, long& req_addr, Request::Type& req_type, unsigned& cpu_id){
    string linebuffer;
   // if(file.eof()){cout << "end of file\n"; return false;}

    if(getline(file, linebuffer)){
		if (linebuffer.length() == 0){
			return false;
		}
        vector<string> split = split_str(linebuffer);
        line l;
        try{
            if(split.size() < 5){
                req_addr = -1;
                return true;
            }
            int type_position = 0;
            for(int i = 0; i < split.size(); i++){
                if((split[i] == "L") || (split[i] == "S") || (split[i] == "P") || (split[i] == "I")){
                    type_position = i;
                    break;
                }
            }
            if((type_position >= 3) && ((type_position+1) < split.size())){
                l.threadID = std::stoi(split[type_position - 3],nullptr,10);
                l.processorID = std::stoi(split[type_position - 2],nullptr,10);
                if(l.processorID > 512){
                    req_addr = -1;
                    return true;
                }
                if(split[type_position - 1] != "-")
                    l.cycleNum = std::stoi(split[type_position - 1],nullptr,10);
                else
                    l.cycleNum = 0;
                l.type = split[type_position];
                l.addr  = std::stoul(split[type_position + 1], nullptr, 10);

                bubble_cnt = l.cycleNum;
                req_addr = l.addr;
                if (l.type =="L") req_type = Request::Type::READ;
                else if (l.type == "S") req_type = Request::Type::WRITE;
                else if (l.type == "I") req_type = Request::Type::INSTRUCTION;

                cpu_id = l.processorID;
            }
            else{
              // i cout << linebuffer << endl;
                req_addr = -1;
                return true;
            }
        }
        catch(std::exception e){
         //   cout << linebuffer << endl;
            req_addr = -1;
            return true;
        }
        return true;
    }
    else{
        return false;
    }	
}



bool Trace::get_rowclone_request(long& bubble_cnt, long& req_addr, Request::Type& req_type)
{
    static bool has_write = false;
    static long write_addr;
    static int line_num = 0;

    static bool batch_cp = false;
    static bool batch_set = false;
    static bool is_src;
    static long src_addr;
    static long src_addr_end;
    static long dst_addr;
    static long dst_addr_end;

    if (has_write){
        bubble_cnt = 0;
        req_addr = write_addr;
        req_type = Request::Type::WRITE;

        has_write = false;
        return true;
    } else if (batch_cp) {
      bubble_cnt = 0;
      if (is_src) {
        req_addr = src_addr;
        req_type = Request::Type::READ;
        if (src_addr < src_addr_end) {
          // FIXME: won't assume cacheline size is 64.
          src_addr += 64;
        }
      } else {
        req_addr = dst_addr;
        req_type = Request::Type::WRITE;
        if (dst_addr < dst_addr_end) {
          // FIXME: won't assume cacheline size is 64.
          dst_addr += 64;
        }
      }
      if (src_addr >= src_addr_end && dst_addr >= dst_addr_end) {
        batch_cp = false;
      } else {
        is_src = !is_src;
      }
      return true;
    } else if (batch_set) {
      bubble_cnt = 0;
      req_addr = dst_addr;
      req_type = Request::Type::WRITE;

      dst_addr += 64;
      if (dst_addr >= dst_addr_end) {
        batch_set = false;
      }
      return true;
    }

    string line;
    getline(file, line);
    line_num ++;
    if (file.eof() || line.size() == 0) {
        file.clear();
        file.seekg(0, file.beg);
        has_write = false;
        line_num = 0;
        return false;
    }

    size_t pos, end;

    char type = line[0];
    switch (type) {
      case 'R': {
        pos = line.find_first_not_of(' ', 1);
        req_addr = stoul(line.substr(pos), &end, 16);
        req_type = Request::Type::READ;

        pos = line.find_first_not_of(' ', pos+end);
        write_addr = stoul(line.substr(pos), &end, 16);
        if (write_addr != 0) {
          has_write = true;
        }

        pos = line.find_first_not_of(' ', pos+end);
        bubble_cnt = stoul(line.substr(pos), NULL, 10);
        break;
      }
      case 'C': {
        batch_cp = true;

        pos = line.find_first_not_of(' ', 1);
        src_addr = stoul(line.substr(pos), &end, 16);

        pos = line.find_first_not_of(' ', pos + end);
        dst_addr = stoul(line.substr(pos), &end, 16);

        pos = line.find_first_not_of(' ', pos + end);
        int cp_data_size = stoul(line.substr(pos), &end, 10);

        pos = line.find_first_not_of(' ', pos + end);
        bubble_cnt = stoul(line.substr(pos), NULL, 10);

        src_addr_end = src_addr + cp_data_size;
        dst_addr_end = dst_addr + cp_data_size;

        req_addr = src_addr;
        is_src = false;
        src_addr += 64;

        req_type = Request::Type::READ;

        break;
      }
      case 'S': {
        batch_set = true;

        pos = line.find_first_not_of(' ', 1);
        dst_addr = stoul(line.substr(pos), &end, 16);

        pos = line.find_first_not_of(' ', pos + end);
        assert(stoul(line.substr(pos), &end, 16) == 0);

        pos = line.find_first_not_of(' ', pos + end);
        int set_data_size = stoul(line.substr(pos), &end, 10);

        pos = line.find_first_not_of(' ', pos + end);
        bubble_cnt = stoul(line.substr(pos), NULL, 10);

        dst_addr_end = dst_addr + set_data_size;

        req_addr = dst_addr;
        req_type = Request::Type::WRITE;
        dst_addr += 64;
        if (dst_addr >= dst_addr_end) {
          batch_set = false;
        }
        break;
      }
    }

    return true;
}
