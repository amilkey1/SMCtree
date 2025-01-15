#pragma once

namespace proj {

    class StopWatch {
    
        public:
            StopWatch() : running(false) {}
            ~StopWatch() {}

            void   start();
            double stop();
            
        private:
            bool    running;
            clock_t started;
            clock_t stopped;
    };

    inline void StopWatch::start() {
        assert(!running);
        started = clock();
        running = true;
    }

    inline double StopWatch::stop() {
        assert(running);
        stopped = clock();
        double seconds = (double)(stopped - started)/CLOCKS_PER_SEC;
        running = false;
        started = stopped = 0L;
        return seconds;
    }

}
