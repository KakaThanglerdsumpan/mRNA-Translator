#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *translate(char *mrna, int num, bool override);
float compare(char *seq1, char *seq2);

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        printf("Command line argument must contain at least one mRNA strand:\n");
        printf("Usage: ./mrna (override) [mRNA1] ([mRNA2] [mRNA3] [...])\n");
        return 1;
    }

    int x = 1;
    bool override = false;

    if (strcmp(argv[1], "override") == 0)
    {
        override = true;
        x = 2;
    }

    if (override == true && argc == 2)
    {
        printf("Command line argument must contain at least one mRNA strand:\n");
        printf("Usage: ./mrna (override) [mRNA1] ([mRNA2] [mRNA3] [...])\n");
        return 1;
    }

    char *mrna[argc - x];
    for (int i = 0; i < argc - x; i++)
    {
        mrna[i] = argv[i + x];
    }

    for (int i = 0; i < argc - x; i++)
    {
        for (int j = 0; j < strlen(mrna[i]); j++)
        {
            if (mrna[i][j] == 'A' || mrna[i][j] == 'U' || mrna[i][j] == 'G' || mrna[i][j] == 'C')
            {
            }
            else
            {
                printf("mRNA strands can only contain A, U, G, and C bases!\nPlease capitalize all of your bases!\n");
                return 1;
            }
        }
    }

    char *aminoAcids[argc - x];
    for (int i = 0; i < argc - x; i++)
    {
        aminoAcids[i] = translate(mrna[i], i + 1, override);
    }

    if (argc - x > 1)
    {
        printf("\nSimilarity:\n");
        for (int i = 0; i < argc - x; i++)
        {
            for (int j = i + 1; j < argc - x; j++)
            {
                float pIdentical = compare(aminoAcids[i], aminoAcids[j]);
                char p = '%';
                printf("\tSeq.%i—Seq.%i: %.2f%c identical\n", i + 1, j + 1, pIdentical, p);
            }
        }
    }

    if (override == false)
    {
        printf("\nTo translate the entire mRNA strand (overriding start and stop codons), please include 'override' as your second command line argument:\n");
        printf("./mrna override [mRNA...]\n");
    }

    for (int i = 0; i < argc - x; i++)
    {
        free(aminoAcids[i]);
    }
    return 0;
}

char *translate(char *mrna, int num, bool override)
{
    int startPos = 0;
    bool foundStart = false;
    bool incCodon = false;

    // array of chars to store amino acid sequence
    char *aminoAcids = malloc((((strlen(mrna) - startPos) / 3) + 1) * sizeof(char));
    memset(aminoAcids, 0, ((strlen(mrna) - startPos) / 3) + 1);

    if (override == false)
    {
        // scans through mrna in windows of three to look for start codon—"AUG"
        for (int i = 0; i <  strlen(mrna); i++)
        {
            char MET[] = "AUG";
            char curWindow[4];
            curWindow[3] = '\0';
            for (int j = i; j < i + 3; j++)
            {
                curWindow[j - i] = mrna[j];
            }
            if (strcmp(curWindow, MET) == 0)
            {
                startPos = i;
                foundStart = true;
                break;
            }
        }
    }

    // once start codon is found, translate each codon into an aminco acid and store it in aminoAcids
    if (foundStart || override == true)
    {
        // look at each codon starting with the start codon
        for (int i = startPos; i < strlen(mrna); i += 3)
        {
            char curCodon[4];
            for (int j = i; j < i + 3; j++)
            {
                if (mrna[j] == '\0')
                {
                    incCodon = true;
                    break;
                }
                curCodon[j - i] = mrna[j];
            }

            // if it is a complete codon, search for the amino acid that it codes for...
            if (incCodon == false)
            {
                if (curCodon[0] == 'A')
                {
                    if (curCodon[1] == 'A')
                    {
                        if ('A' == curCodon[2] || 'G' == curCodon[2])
                        {
                            aminoAcids[(i - startPos) / 3] = 'K';
                        }
                        else
                        {
                            aminoAcids[(i - startPos) / 3] = 'N';
                        }
                    }
                    if (curCodon[1] == 'U')
                    {
                        if (curCodon[2] == 'G')
                        {
                            aminoAcids[(i - startPos) / 3] = 'M';
                        }
                        else
                        {
                            aminoAcids[(i - startPos) / 3] = 'I';
                        }
                    }
                    if (curCodon[1] == 'C')
                    {
                        aminoAcids[(i - startPos) / 3] = 'T';
                    }
                    if (curCodon[1] == 'G')
                    {
                        if ('A' == curCodon[2] || 'G' == curCodon[2])
                        {
                            aminoAcids[(i - startPos) / 3] = 'R';
                        }
                        else
                        {
                            aminoAcids[(i - startPos) / 3] = 'S';
                        }
                    }
                }
                if (curCodon[0] == 'U')
                {
                    if (curCodon[1] == 'A')
                    {
                        if ('A' == curCodon[2] || 'G' == curCodon[2])
                        {
                            aminoAcids[(i - startPos) / 3] = '*';
                            if (override == false)
                            {
                                break;
                            }
                        }
                        else
                        {
                            aminoAcids[(i - startPos) / 3] = 'Y';
                        }

                    }
                    if (curCodon[1] == 'U')
                    {
                        if ('A' == curCodon[2] || 'G' == curCodon[2])
                        {
                            aminoAcids[(i - startPos) / 3] = 'L';
                        }
                        else
                        {
                            aminoAcids[(i - startPos) / 3] = 'F';
                        }
                    }
                    if (curCodon[1] == 'C')
                    {
                        aminoAcids[(i - startPos) / 3] = 'S';
                    }
                    if (curCodon[1] == 'G')
                    {
                        if (curCodon[2] == 'A')
                        {
                            aminoAcids[(i - startPos) / 3] = '*';
                            if (override == false)
                            {
                                break;
                            }
                        }
                        else if (curCodon[2] == 'G')
                        {
                            aminoAcids[(i - startPos) / 3] = 'W';
                        }
                        else
                        {
                            aminoAcids[(i - startPos) / 3] = 'C';
                        }
                    }
                }
                if (curCodon[0] == 'C')
                {
                    if (curCodon[1] == 'A')
                    {
                        if ('A' == curCodon[2] || 'G' == curCodon[2])
                        {
                            aminoAcids[(i - startPos) / 3] = 'Q';
                        }
                        else
                        {
                            aminoAcids[(i - startPos) / 3] = 'H';
                        }
                    }
                    if (curCodon[1] == 'U')
                    {
                        aminoAcids[(i - startPos) / 3] = 'L';
                    }
                    if (curCodon[1] == 'C')
                    {
                        aminoAcids[(i - startPos) / 3] = 'P';
                    }
                    if (curCodon[1] == 'G')
                    {
                        aminoAcids[(i - startPos) / 3] = 'R';
                    }
                }
                if (curCodon[0] == 'G')
                {
                    if (curCodon[1] == 'A')
                    {
                        if ('A' == curCodon[2] || 'G' == curCodon[2])
                        {
                            aminoAcids[(i - startPos) / 3] = 'E';
                        }
                        else
                        {
                            aminoAcids[(i - startPos) / 3] = 'D';
                        }
                    }
                    if (curCodon[1] == 'U')
                    {
                        aminoAcids[(i - startPos) / 3] = 'V';
                    }
                    if (curCodon[1] == 'C')
                    {
                        aminoAcids[(i - startPos) / 3] = 'A';
                    }
                    if (curCodon[1] == 'G')
                    {
                        aminoAcids[(i - startPos) / 3] = 'G';
                    }
                }
            }
        }
        printf("Amino Acid Sequence %i: ", num);
        for (int i = 0; i < (((strlen(mrna) - startPos) / 3) + 1); i++)
        {
            if (aminoAcids[i] != '\0')
            {
                printf("%c ", aminoAcids[i]);
            }
        }
        if (incCodon == true)
        {
            printf("———INCOMPLETE CODON———");
        }
        printf("\n");
    }
    else
    {
        printf("Amino Acid Sequence %i: NO START CODON FOUND!\n", num);
    }
    return aminoAcids;
}

float compare(char *seq1, char *seq2)
{
    int length1 = strlen(seq1);
    int length2 = strlen(seq2);
    float divisor; //divisor for finding percentage of similarity
    float nIdentical = 0; // number of identical amino acids
    float percent; // percent similarity between the two amino acid sequences
    int limit;

    if (length1 == 0 && length2 == 0)
    {
        return 0;
    }
    if (length1 > length2)
    {
        divisor = length1;
        limit = length2;
    }
    else if (length1 < length2)
    {
        divisor = length2;
        limit = length1;
    }
    else
    {
        divisor = length1;
        limit = length1;
    }

    for (int i = 0; i < limit; i++)
    {
        if (seq1[i] == seq2[i])
        {
            nIdentical++;
        }
    }
    percent = (nIdentical / divisor) * 100;
    return percent;
}
