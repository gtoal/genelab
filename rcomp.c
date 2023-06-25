#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv)
{
  int c;
  char *s, *dna;
  char *t;
  if (argc != 2) {
    fprintf(stderr, "syntax: rcomp GATTACA\n");
    exit(1);
  }
  dna = argv[1];
  t = strdup(dna);*t++ = '\0';
  s = dna;
  while (*s != '\0') { // complement
    c = *s++;
    if (c == 'A' || c == 'a') {
      c = 'T';
    } else if (c == 'T' || c == 't') {
      c = 'A';
    } else if (c == 'G' || c == 'g') {
      c = 'C';
    } else if (c == 'C' || c == 'c') {
      c = 'G';
    } else {
      fprintf(stderr, "rcomp: expected one of [ACGT] at %s\n", --s);
      exit(1);
    }
    *t++ = c;
  }
  for (;;) { // reverse
    c = *--t;
    if (c == '\0') break;
    putchar(c);
  }
  putchar('\n');
  exit(0);
  return 0;
}
