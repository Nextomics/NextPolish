#include "scorechain.h"

PolishResult * score_chain(const char * tigname, Configure* configure)
{
	Contig* contig = contig_init(tigname, configure, sizeof(Kmer));

	contig->read_fliter = contig_read_fliter1;
	contig_score_correct(contig, 0, contig->length - 1, 0x1, contig->configure->indel_balance_factor_sgs);

	PolishResult* result = contig_get_contig(contig, 0, contig->length - 1, FLAG_ZERO | FLAG_COVERAGE);

	contig_destory(contig);

	return result;
}

PolishResult * td_score_chain1(const char * tigname, Configure* configure)
{
	Contig* contig = contig_init(tigname, configure, sizeof(Kmer));

	contig->read_fliter = contig_read_fliter2;
	contig_score_correct(contig, 0, contig->length - 1, 0x1, contig->configure->indel_balance_factor_lgs);

	PolishResult* result = contig_get_contig(contig, 0, contig->length - 1, 0);

	contig_destory(contig);

	return result;
}
