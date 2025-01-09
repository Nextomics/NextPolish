#ifndef SCORE_CHAIN_H
#define SCORE_CHAIN_H

#include "contig.h"

PolishResult* score_chain(const char* tigname, Configure* configure);

PolishResult* td_score_chain1(const char* tigname, Configure* configure);

#endif
