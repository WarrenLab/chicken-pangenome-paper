#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN_COVERAGE 10
#define MAX_COVERAGE 100
#define MAX_SEQNAME_LENGTH 200
#define MAX_PILEUP_SIZE 8001
#define min(a, b) a < b ? a : b

float calcRefFreq(char *pileup) {
  float minor_allele_count;
  char c;

  int ref_bases = 0;
  int nonref_A = 0;
  int nonref_C = 0;
  int nonref_G = 0;
  int nonref_T = 0;
  int i = 0;

  for (c = pileup[0]; c != '\0'; c = pileup[++i]) {
    if (c == '.' || c == ',') {
      ref_bases++;
    }

    else {
      switch (toupper(c)) {
      case 'A':
        nonref_A++;
        break;
      case 'C':
        nonref_C++;
        break;
      case 'G':
        nonref_G++;
        break;
      case 'T':
        nonref_T++;
        break;
      }
    }
  }

  int total_bases = ref_bases + nonref_A + nonref_C + nonref_G + nonref_T;
  int max_nonref = 0;

  if (nonref_A > max_nonref) {
    max_nonref = nonref_A;
  }
  if (nonref_C > max_nonref) {
    max_nonref = nonref_C;
  }
  if (nonref_G > max_nonref) {
    max_nonref = nonref_G;
  }
  if (nonref_T > max_nonref) {
    max_nonref = nonref_T;
  }

  minor_allele_count = min(max_nonref, ref_bases);
  if ((max_nonref + ref_bases) * 1.0 / total_bases > 0.9 &&
      minor_allele_count / (ref_bases + max_nonref) > 0.25) {
    return max_nonref * 1.0 / (ref_bases + max_nonref);
  }

  return -1.0;
}

int main(int argc, char *argv[]) {
  int position, coverage;
  char base;
  float refFreq;

  char *seqName = malloc(MAX_SEQNAME_LENGTH * sizeof(char));
  char *pileup = malloc(MAX_PILEUP_SIZE * sizeof(char));
  char *quals = malloc(MAX_PILEUP_SIZE * sizeof(char));

  while (scanf("%s\t%d\t%c\t%d\t%s\t%s", seqName, &position, &base, &coverage,
               pileup, quals)) {
    if (coverage >= MIN_COVERAGE && coverage <= MAX_COVERAGE) {
      refFreq = calcRefFreq(pileup);
      if (refFreq > 0.0) {
        printf("%s\t%d\t%c\t%d\t%s\t%f\n", seqName, position, base, coverage,
               pileup, refFreq);
      }
    }
  }
  return 0;
}
