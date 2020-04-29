#ifndef PROGRESSBAR_HPP
#define PROGRESSBAR_HPP

#include <chrono>
#include <stdio.h>

class ProgressBar {
private:
    unsigned int ticks = 0;
    const unsigned int total_ticks;
    const unsigned int bar_width;
    const char complete_char = '=';
    const char incomplete_char = ' ';
    const std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

public:
    ProgressBar(unsigned int total, unsigned int width, char complete, char incomplete) :
            total_ticks {total},
            bar_width {width},
            complete_char {complete},
            incomplete_char {incomplete} {}

    ProgressBar(unsigned int total, unsigned int width) :
            total_ticks {total},
            bar_width {width} {}

    unsigned int operator++() { return ++ticks; }

    void display() const
    {
        using std::chrono::steady_clock;
        using std::chrono::duration_cast;
        using std::chrono::milliseconds;

        double progress = (double)ticks/total_ticks;
        unsigned int pos = (unsigned int)(bar_width*progress);

        steady_clock::time_point now = steady_clock::now();
        double time_elapsed = (double)(duration_cast<milliseconds>(now-start_time).count())/1000.0;
        double time_remain = time_elapsed/progress*(1 - progress);

        // output
        fprintf(stdout, "\r%3d%%[", int(progress*100.0));
        for (unsigned int i = 0; i < bar_width; ++i) {
            if (i < pos) putc(complete_char, stdout);
            else if (i == pos) putc('>', stdout);
            else putc(incomplete_char, stdout);
        }
        fprintf(stdout, "] %.3fs/%.3fs", time_elapsed, time_elapsed+time_remain);
        fflush(stdout);
    }

    void done() const
    {
        display();
        putc('\n', stdout);
    }
};

#endif // PROGRESSBAR_HPP
