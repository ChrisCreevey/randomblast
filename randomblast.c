#define STD_CODON_NUM 300			/* standard length of the gene in codons */
#define maxnamlen 50			/* maximum length of the gene name */


#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

struct sequence{

	char *bases;					/* The sequence which is stored in dynamically allocted memory  */
	char name[maxnamlen + 1 ];	/* the name of the sequence */
	int seq_num;				/* the number of the sequence in the file */
	int tag;					/* A tag for the seguence */
	struct sequence *next;		/* pointer to the next sequence */
	struct sequence *previous;	/* pointer to the previous sequence */

	} list_entry;



int print_seq(int seq);
int delete_sequence(char *name);
void totext(int c, char *array);
char inttotext(int c);
void open_input_file(void);
int read_file (int fastaformat);
char read_sequence(int seq_num, struct sequence *new);
void clean_exit(void);
void clear_memory(void);




struct sequence *start ='\0';  /*pointer to first sequence in the list */
struct sequence *last = '\0'; /*point to last sequece in the list */

FILE *file = '\0', *summary = '\0', *fastafile = '\0';
int num_of_seqs = 0, blastnum = 0;
char filename[30], **todelete = '\0', blast_query[maxnamlen];



int main(void)
	{
	struct sequence *position = start;
	char array[30], string[100], c = '\0', tmp[100];
	int i = 0, deleted = 0, candidate = 0, j=0,blast_attempts = 0, system_result = 1;
	
	/*seed the rand number with the calander time ---- */
	  srand((unsigned) time(NULL)/2);
	blast_query[0] = '\0';
	tmp[0] = '\0';
	summary = fopen("summary.txt", "w");
	system("date >> date.txt");
	open_input_file();
        printf("number of sequences = %d\n", num_of_seqs);
	todelete = malloc(num_of_seqs * sizeof(char *));
	if(todelete == '\0')
	  { 
	  printf("out of memory\n");
	  clean_exit();
	  }
	for(i=0; i<num_of_seqs; i++) 
	  {
	    todelete[i] = malloc(100 *sizeof(char));
	    if(todelete[i] == '\0')
	      {
		printf("out of memory\n");
		clean_exit();
	      }
	  }
	 /* system call to formatdb to prepare the database for searching  
	 THIS IS DONE ONLY ONCE */
		
	/* if(system("formatdb -i dbaa -p") != 0)
 		{
		printf("Error in calling formatdb\n");
		clean_exit();
	    }
	  */  
	    
	while(num_of_seqs > 0)
		{
		if(blastnum != 0)
			{
			totext(blastnum-1, array);
			strcpy(filename, array);
			strcat(filename, ".blast");
			fprintf(summary, "%s\t",filename); 
			if((file = fopen(filename, "r")) != '\0')		/* check to see if the file is there */
				{

				/* if it is */
				/* open an output file that will contain the sequences from the blast result file */
				strcpy(filename, array);
				strcat(filename, "fas.txt");
				fastafile = fopen(filename, "w");

				/* read down to line 17 of the result file to see if there are any results */
				  i = 0;j=0;
				while(i < 16 && !feof(file))
					{
					c = getc(file);
					if(c == '\r' || c == '\n')
					    i++;
					j++;
					}
				fscanf(file, "%s", string);
			       
				/* check to see if this line has the phrase "no hits were found"  */
			       
				if(strncmp(string,"*****", 100) != 0)   /* if there were hits */
					{
					for(j=0; j<num_of_seqs; j++) todelete[j][0] = '\0';
					
					while(i< 19 && !feof(file))
						{
						c = getc(file);
						if(c == '\r' || c == '\n') i++;
						}
						
					fscanf(file, "%s", string);
					strncpy(string, string, maxnamlen);
				       
					deleted = 0; c = '\0';
					while(!feof(file) && c != '\r' && c != '\n')
						{
                                            
						strcpy(todelete[deleted], string);
					        deleted++;
						while(!feof(file) && (c = getc(file))!= '\n' && c != '\r');
						if((c = getc(file)) != '\r' && c != '\n')
						  {
						    i =0;
						    while(c != ' ')
						      {
								if(i<maxnamlen)
									{
									string[i] = c;
									i++;
									}
							  c = getc(file);
						      }
						    string[i] = '\0';
						    
						  
						  }
						}
					fclose(file);
			/*		if(deleted >= 0)
					  {
					    strcpy(tmp, "rm ");
					    strcat(tmp, filename);
					   
					    system(tmp); 
					  }
			*/			
					/* next step is to delete these sequences from memory */
					position = start; i = 0;
					while(position != '\0' && position->tag == FALSE)
						position = position->next;
					
					while(position != '\0' && i != candidate) /*** start by making sure that the candidate is printed out (BLAST might not always find the query sequence) **/
						{
					    position = position->next;
					    if(position->tag == TRUE) i++;
						}
				       
					untag_sequence(position->name);
				            
				       
					for(i=0; i<deleted; i++)
						{
						if(strncmp(blast_query, todelete[i], maxnamlen) != 0)  /* make sure that this sequence is not the same as the query (because we have already printed that out) */
							untag_sequence(todelete[i]);
						}

					}
				    else
				        { /* if there were no hits then we need to delete the sequence that we used to blast the last time..... candiate */   
					deleted = 0;
					position = start; i = 0;
				        while(position != '\0' && i != candidate)
				           {
					    position = position->next;
					    if(position->tag == TRUE) i++;
				           }
				       
				       /* if(delete_sequence(position->name) == TRUE)  */
				       
				       	untag_sequence(position->name);
  
				        }
				fprintf(summary, "%d\n", deleted);
				fclose(fastafile);
				
				}
			        else
				  printf("Couldn't open file %s\n",filename);
				
			}
		if(num_of_seqs > 0)
		  {
		    /* randomly pick a sequence file from memory */
		    candidate = (int)fmod(rand(), num_of_seqs);
		    if(!print_seq(candidate))
			{
			printf("WARNING.....could not find next sequence to blast\n");
			clean_exit();
			}
		
		    /* Write everyhing in memory to the data base file */
			if((fmod((float)blastnum, 500.0)) == 0)
				{
				file = fopen("database.txt", "w");
				position = start;
				while(position != '\0')
					{
					if(position->tag== TRUE)
						{
						fprintf(file, ">%s\n", position->name);
						fprintf(file, "%s\n", position->bases);
						}
					position = position->next;
					}
				fflush(file);
				fclose(file);
				
			
				/* system call to formatdb to prepare the database for searching */
			
				if(system("formatdb -i database.txt -p") != 0)
				  {
				printf("Error in calling formatdb\n");
				clean_exit();
				  }
				}
					
		
		    /* system call to blast the input file against the database */
		    strcpy(string, "blastall -p blastp -d database.txt -i blast.in -e 0.0000001 -v 10000 -b 1 -o ");
		    totext(blastnum, array);
		    strcat(string, array);
		    strcat(string, ".blast");
                    blast_attempts = 0;
                    do
                        {
                        system_result = system(string);
                        blast_attempts++;
                        }while(system_result != 0 && blast_attempts < 100);  /* The program will attempt to run the blast 100 times before quitting */
		    if(system_result !=0)
			{
			printf("Error in calling blastall\n");
                        /* print the remaining sequences to the file database.txt*/
                        file = fopen("database.txt", "w");
                        position = start;
                        while(position != '\0')
                            {
                            if(position->tag== TRUE)
                                {
                                fprintf(file, ">%s\n", position->name);
                                fprintf(file, "%s\n", position->bases);
                                }
                            position = position->next;
                            }
                        fflush(file);
                        fclose(file);
			clean_exit();
			}
		  }	
		blastnum++;
		printf("Number of seqs = %d, blast number = %d\n", num_of_seqs, blastnum);
		}
	fclose(summary);
	system("date >> date.txt");
	clear_memory();
	printf("Finished!\n");
	}
	


int print_seq(int seq)
	{
	struct sequence * position = start;
	int i=0, found = FALSE;
	
	while(position != '\0' && position->tag == FALSE)
		position = position->next;
	
	while(position != '\0' && i != seq)
		{
		position = position->next;
		if(position->tag == TRUE) i++;
		}
	if(i == seq && position != '\0')
		{
		found = TRUE;
		file = fopen("blast.in", "w");
		strncpy(blast_query, position->name, maxnamlen);
		fprintf(file, ">%s\n", position->name);
		fprintf(file, "%s\n", position->bases);
		fclose(file);
		}
	return(found);
	}
	




int delete_sequence(char *name)
	{
	struct sequence * position = start;
	int found = FALSE;
	
	while(position != '\0')
		{
		if(strncmp(name, position->name, 100) == 0) /*this is the sequence to delete */
			{
			/* first print the sequence to the fasta file for this blast result*/
			fprintf(fastafile, ">%s\n", position->name);
			fprintf(fastafile, "%s\n", position->bases);
			/* now delete the sequence from memory */
			found = TRUE;
			if(position->previous != '\0') (position->previous)->next = position->next;
			if(position->next != '\0') (position->next)->previous = position->previous;
			if(position == start) start = position->next;
			if(position == last) last = position->previous;
			free(position->bases);
			free(position);
			position = '\0';
			
			}
		if(position != '\0')position = position->next;
		}
	if(!found) printf("Cannot find sequence name %s in memory\n", name);
	return(found);
	}

int untag_sequence(char *name)
	{
	struct sequence * position = start;
	int found = FALSE, present = FALSE;
	
	while(position != '\0' && !found && !present)
		{
		if(strncmp(name, position->name, maxnamlen) == 0) /*this is the sequence to delete */
			{
                        present = TRUE;
			/* first print the sequence to the fasta file for this blast result*/
			if(position->tag != FALSE)
				{
				fprintf(fastafile, ">%s\n", position->name);
				fprintf(fastafile, "%s\n", position->bases);
				/* now untag the sequence in memory */
				found = TRUE;
				position->tag = FALSE;
				num_of_seqs--;
				}
			}
		if(position != '\0' && !present) position = position->next;
		}
        if(!present) printf("ERROR: cannot find sequence named %s\n", name);

	return(found);
	}
	




void totext(int c, char *array)
	{
	int count = 0, i = 0;
	char tmp[30];
	
	while(pow(10,count) <= c)
		{
		tmp[count] = inttotext(fmod(c/pow(10, count), 10));
		count++;
		}
	if(c == 0)
	  {
	    count++;
	    tmp[0] = '0';
	  }
	tmp[count] = '\0';
	array[count] = '\0';
	for(i= count-1; i>=0; i--)
		{
		array[(count-1) - i] = tmp[i];
		}
	}
	
	
char inttotext(int c)
	{
	
	switch(c)
		{
		case 1 :
			return('1');
			break;
		case 2 :
			return('2');
			break;
		case 3 :
			return('3');
			break;
		case 4 :
			return('4');
			break;
		case 5 :
			return('5');
			break;
		case 6 :
			return('6');
			break;
		case 7 :
			return('7');
			break;
		case 8 :
			return('8');
			break;
		case 9 :
			return('9');
			break;
		default:
			return('0');
			break;
		}
		
	}
		



void open_input_file(void)

	{
	int fastaformat = TRUE;


	if((file = fopen("dbaa", "r")) == '\0')		/* check to see if the file is there */
		{
		printf("\n\n\tCannot open the sequence file, named %s\n\n", filename);
		}

	printf("opened %s.....\n", filename);
			
	fastaformat = read_file(fastaformat);
	if(!fastaformat)
		{
		printf("%s is not in fasta format\n");
		clean_exit();
		}
	fclose(file);

	}










int read_file (int fastaformat)
	{
	
	struct sequence *new = '\0';
	int j;
	char c;
	
	if((c = getc(file)) == '>')
		{
		fastaformat = TRUE;
		j = 0;
		do{
				new = malloc(sizeof(list_entry));
				if(!new)
					{
					printf("\n\t Out of memory\n");
					clean_exit();
					}
			
			c = read_sequence(j, new);
			j++;
			
			}while((c == '>') && (feof(file) != 1));
		}
			
	else
		{
		/* if the first character in the file is not a '>' then  we assume the format is wrong */
 		printf("\n\n\tThis does not seem   to be a fasta formatted file\n\n");
		fastaformat = FALSE;
		}
	num_of_seqs = j;	
	return(fastaformat);
	
	}
			
	
	
	
	
	
	
char read_sequence(int seq_num, struct sequence *new)
	{
	
	char c = '\0';
	int i = 0, j = 0, place = 0, value = 0, memory_allocations = 0;
	


	/* assign the sequence number */
	new->seq_num = seq_num;
	

	/* initialise the tag (equal to TRUE) */
	new->tag = TRUE;
	

	/* read in the name of the sequence */
	j = 0;
	c = getc(file);
	while(c == ' ') c = getc(file);
	do{				
	 
	 	new->name[j] = c;
		if(j < maxnamlen) j++;
		

		}while((c = getc(file)) != '\n' && c != '\r' && c != ' ' && feof(file) == 0);  /*read up until the first space */
	new->name[j] = '\0';				/* append a '\0' terminator */
        if(c == ' ')
            {
            while( c != '\n' && c != '\r' && !feof(file))
                c = getc(file);
            }
            
	/* read in the sequence */
    i = 0;
    j = 0;
    
    new->bases = malloc((STD_CODON_NUM * (sizeof(char))));   /* Allocate the memory to store the sequence in codons */
	memory_allocations = 1;
	c = getc(file);


	do{			
		if(c != '\n' && c != '\r' && feof(file) == 0 && c != ' ')  /* if not an end of line or end of file */
			{
			new->bases[i] = c;
			i++;			
			}
		
		if(i >= ((STD_CODON_NUM * memory_allocations) - 1))  /* if the sequence has exceeded the length of STD_CODON_NUM * the number of memory allocations */
			{
			memory_allocations++;
			new->bases = realloc(new->bases, ((STD_CODON_NUM * (sizeof(char))) * memory_allocations));    /* reallocate the memory space to the size required to fit the STD_CODON_NUM * new memory allocations */
			if(!new->bases)
				{
				printf("\n\t Out of memory\n");
				clean_exit();
				}
			}

		}while(((c= getc(file)) != '>') && (feof(file) == 0));

		
	new->bases[i] = '\0';	/* NULL is the terminator value for the sequence */
	

	
	/* now assign the pointers to the previous and next sequence */
	if (!last || !start) 	/* If this is a new list */
		{
		last = new;
		start = new;
		}		
	else (last)->next = new;
	new->next = '\0';
	new->previous = last;
	last = new;

	return(c);


	}

void clean_exit(void)
	{
	clear_memory();
	exit(1);
	}




void clear_memory(void)
	{
	int i;
	struct sequence *new = '\0';

	if(start)
		{
		while(start) 
			{
			 new = start->next;
			 free(start->bases);
			 free(start); 
			 start = new;
			}
		start = '\0';
		last = '\0';
		}
	}


