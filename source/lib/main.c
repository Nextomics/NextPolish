#include "scorechain.h"
#include "kmercount.h"
#include "snpphase.h"
#include "snpvalid.h"
#include "lgspolish.h"

const char* argvs[] = { "scorechain", "kmercount", "snpphase", "snpvalid", "lgspolish" };
contig_run func[] = { score_chain,kmer_count,snp_phase,snp_valid,lgspolish};

int32_t check_argv(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	int32_t step = check_argv(argc, argv);
	if (step) {
		if (step == 5) {
			Configure* configure = config_init(argv[2], NULL, argv[3]);
			contig_total(func[step - 1], step, configure, -1);
		}
		else {
			Configure* configure = config_init(argv[2], argv[3], argv[4]);
			contig_total(func[step - 1], step, configure, -1);
		}
	}
	return 0;
}


int32_t check_argv(int argc, char *argv[]) {
	int32_t step = -1;
	char* command =
		"Usage: %s <command> [options]\n\n\
Commands:\n\
	scorechain		score chain run\n\
				eg. scorechain fastafn sgsbamf > output.fa\n\
	kmercount		kmer count run\n\
				eg. kmercount fastafn sgsbamf > output.fa\n\
	snpphase		snp phase run\n\
				eg. snpphase fastafn sgsbamf lgsbamf > output.fa\n\
	snpvalid		snp valid run\n\
				eg. snpvalid fastafn sgsbamf > output.fa\n\
	lgspolish 		lgs polish run\n\
				eg. lgspolish fasta_file lgsbamf > output.fa\n\n";

	if (argc > 1) {
		for (int i = 0; i < sizeof(argvs) / sizeof(const char *); i++) {
			if (strcmp(argvs[i], argv[1]) == 0) {
				step = i;
				break;
			}
		}
	}
	step++;
	switch (step) {
	case 1:
	case 2:
	case 4:
	case 5:
		if (argc != 4) {
			printf("%s %s fastafn lgsbam\n", argv[0], argv[1]);
			step = 0;
		}
		break;
	case 3:
		if (argc != 5) {
			printf("%s %s fastafn bamfn thirdbamfn\n", argv[0], argv[1]);
			step = 0;
		}
		break;
	default:
		printf(command, argv[0]);
		step = 0;
	}
	return step;
}
