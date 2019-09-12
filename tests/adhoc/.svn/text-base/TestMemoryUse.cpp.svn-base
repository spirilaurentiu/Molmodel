#include "Molmodel.h"
#include <iostream>

#if defined (__linux__)
# include <sys/sysinfo.h>

#elif defined(__APPLE__)
# include <mach/task.h>
# include <mach/mach_init.h>

#elif defined(_WINDOWS)
# include <windows.h>
# include "psapi.h"

#else
# include <sys/resource.h>
#endif

/// The amount of memory currently being used by this process, in bytes.
/// By default, returns the full virtual arena, but if resident=true,
/// it will report just the resident set in RAM (if supported on that OS).
size_t memory_used (bool resident=false)
{
#if defined(__linux__)
    // Ugh, getrusage doesn't work well on Linux.  Try grabbing info
    // directly from the /proc pseudo-filesystem.  Reading from
    // /proc/self/statm gives info on your own process, as one line of
    // numbers that are: virtual mem program size, resident set size,
    // shared pages, text/code, data/stack, library, dirty pages.  The
    // mem sizes should all be multiplied by the page size.
    size_t size = 0;
    FILE *file = fopen("/proc/self/statm", "r");
    if (file) {
        unsigned long vm = 0;
        fscanf (file, "%ul", &vm);  // Just need the first num: vm size
        fclose (file);
       size = (size_t)vm * getpagesize();
    }
    return size;

#elif defined(__APPLE__)
    // Inspired by:
    // http://miknight.blogspot.com/2005/11/resident-set-size-in-mac-os-x.html
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    task_info(current_task(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
    size_t size = (resident ? t_info.resident_size : t_info.virtual_size);
    return size;

#elif defined(_WINDOWS)
    // According to MSDN...
    PROCESS_MEMORY_COUNTERS count;
    if (GetProcessMemoryInfo (GetCurrentProcess(), &count, sizeof (count)))
        return count.PagefileUsage;
    else return 0;

#else
    // No idea what platform this is
    return 0;   // Punt
#endif
}

#include <cstring>
#include <ctime>

using namespace SimTK;
using namespace std;

void testRnaResources(int sequenceLength)
{
    time_t initialTime; time(&initialTime);
    size_t initialMemFootprint = memory_used();

    char resBuf[100];
    RNA rna("", false);
    for (int b = 0; b < sequenceLength; ++b) {
        // itoa(b+1, resBuf, 10);
        sprintf(resBuf, "%d", b+1);
        rna.appendResidue(resBuf, RibonucleotideResidue::Adenylate().withPhosphodiester());
    }

    time_t compoundTime; time(&compoundTime);
    size_t compoundMemFootprint = memory_used();

    CompoundSystem system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem forces(system);
    DuMMForceFieldSubsystem dumm(system);
    dumm.loadAmber99Parameters();
    dumm.setAllGlobalScaleFactors(0);
    dumm.setBondTorsionGlobalScaleFactor(1);
    rna.assignBiotypes();
    system.adoptCompound(rna);

    time_t preModelTime; time(&preModelTime);
    system.modelCompounds();
    time_t postModelTime; time(&postModelTime);

    State state = system.realizeTopology();
    system.realize(state, Stage::Position);
    VerletIntegrator integrator(system, 1e-3);
    TimeStepper timeStepper(system, integrator);
    timeStepper.initialize(state);
    timeStepper.stepTo(0.020);

    time_t simTime; time(&simTime);
    size_t systemMemFootprint = memory_used();

    printf("Compound: %.2f Mb/%d bases", (compoundMemFootprint - initialMemFootprint)/1e6, sequenceLength);
    printf("\tSystem: %.2f Mb/%d bases\n", (systemMemFootprint - compoundMemFootprint)/1e6, sequenceLength);

    printf("Compound: %.2f s/%d bases", difftime(preModelTime, initialTime), sequenceLength);
    printf("\tModel: %.2f s/%d bases", difftime(postModelTime, preModelTime), sequenceLength);
    printf("\tSimulation: %.2f s/%d bases /0.020 ps\n", difftime(simTime, postModelTime), sequenceLength);
}

void testMemoryUse() 
{
    // testRnaResources(1);
    testRnaResources(10);
    // testRnaResources(100);
    // testRnaResources(200);
    testRnaResources(400);
    // testRnaResources(800);
    // testRnaResources(1600);
    // testRnaResources(2000);
    // testRnaResources(2200);
    // testRnaResources(2300);
    // testRnaResources(3000);
    // testRnaResources(6000);
    // testRnaResources(6200);
    // testRnaResources(6250);
    // testRnaResources(6260);
    // testRnaResources(6270);
    // testRnaResources(6280);
    // testRnaResources(6282);
    // testRnaResources(6284);
    // testRnaResources(6285);
    // testRnaResources(6286);
    // testRnaResources(6288);
    // testRnaResources(6290);
    // testRnaResources(6300);
    // testRnaResources(6350);
    // testRnaResources(6400);
    // testRnaResources(6450);
    // testRnaResources(6500);
    // testRnaResources(6700);
    // testRnaResources(7000);
    // testRnaResources(8000);
    // testRnaResources(9000);
    // testRnaResources(10000);
    // testRnaResources(50000);
}

int main() {
    try {
        testMemoryUse();
        cout << "PASSED" << endl;
        return 0;
    }
    catch (const std::exception& e)
    {
        printf("EXCEPTION THROWN: %s\n", e.what());
        return 1;
    }
    catch (...)
    {
        printf("UNKNOWN EXCEPTION THROWN\n");
    }    return 1; 
}
