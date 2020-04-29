#ifndef _IO_H_
#define _IO_H_

#include <string>
#include <stdio.h>

#include "ini.h"

#ifndef WIN32
#define ANSI_RESET        "\x1b[0m"
#define ANSI_BOLD         "\x1b[1m"
#define ANSI_BOLD_RED     "\x1b[1;31m"
#define ANSI_BOLD_YELLOW  "\x1b[1;33m"
#define ANSI_BOLD_BLUE    "\x1b[1;34m"
#else
#define ANSI_RESET
#define ANSI_BOLD
#define ANSI_BOLD_RED
#define ANSI_BOLD_YELLOW
#define ANSI_BOLD_BLUE
#endif

#define MSG_WAR ANSI_BOLD_YELLOW "warning: " ANSI_RESET
#define MSG_ERR ANSI_BOLD_RED "error: " ANSI_RESET

void skip_nline(FILE *fp, int n);

// skip header of RINEX format files
void skip_header(FILE *fp);

//void skip_nline(std::ifstream &in, int n);

void backspace(FILE *fp);

#endif
