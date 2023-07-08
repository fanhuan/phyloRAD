#include <iostream>
#include <string>
#include <gzstream.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <fstream>
#include "stdlib.h"
#include "memory.h"
using namespace std;
const size_t Max_Kmer_Size=100;

/*
the latest revision:
There is a -MinMatch parameter, if the value is >=2, each end is required to have one match,
or there are at least this number of non-consecutive matches in a single read.
*/


struct kmer_t
{
	uint64_t kmer;
};


struct bucket
{
	struct kmer_t kmer_t;	//32 bp
	bucket *nxt_bucket;
};


struct hashtable
{
	struct bucket **store_pos;
	size_t ht_sz;
};

struct read_t
{
	char tag[1000];
	bool error_nt[1000];
	char c_seq[1000];//char representation

	uint64_t read_bits[100];//bit representation
	//char read[1000];//char representation
	int readLen;// read length
};

static inline uint64_t get_rev_comp_seq(uint64_t seq, int seq_size)
{
	seq =~seq;

	seq = ((seq & 0x3333333333333333 )<< 2) | ((seq & 0xCCCCCCCCCCCCCCCC )>> 2);
	seq = ((seq & 0x0F0F0F0F0F0F0F0F )<< 4) | ((seq & 0xF0F0F0F0F0F0F0F0 )>> 4);
	seq = ((seq & 0x00FF00FF00FF00FF )<< 8) | ((seq & 0xFF00FF00FF00FF00 )>> 8);
	seq = ((seq & 0x0000FFFF0000FFFF )<<16) | ((seq & 0xFFFF0000FFFF0000 )>>16);
	seq = ((seq & 0x00000000FFFFFFFF )<<32) | ((seq & 0xFFFFFFFF00000000 )>>32);

	return seq >> (64 - (seq_size*2));
}


uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed )
{
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);

	while(data != end)
	{
		uint64_t k = *data++;

		k *= m;
		k ^= k >> r;
		k *= m;

		h ^= k;
		h *= m;
	}

	const unsigned char * data2 = (const unsigned char*)data;

	switch(len & 7)
	{
	case 7: h ^= uint64_t(data2[6]) << 48;
	case 6: h ^= uint64_t(data2[5]) << 40;
	case 5: h ^= uint64_t(data2[4]) << 32;
	case 4: h ^= uint64_t(data2[3]) << 24;
	case 3: h ^= uint64_t(data2[2]) << 16;
	case 2: h ^= uint64_t(data2[1]) << 8;
	case 1: h ^= uint64_t(data2[0]);
	        h *= m;
	};

	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
}


 uint64_t* str2bitsarr(const char * c_str,int len, uint64_t* b_str,int arr_sz )
{

	for (int  k=0;k<arr_sz;++k)
	{
		b_str[k]=0;
	}
	int arr_sz_needed=len/32+1;
	int rem=len%32;
	if(rem==0)
	{arr_sz_needed--;}

	int beg_arr_idx=arr_sz-arr_sz_needed;
	if(rem==0&&arr_sz_needed>0)
	{rem=32;}
	for (int k=0;k<len;k++)
	{
		if(rem==0)
		{beg_arr_idx++;rem=32;}


		switch(c_str[k])
		{
		case ('A'):case ('a'):case ('0'):
			b_str[beg_arr_idx]<<=2;//L_shift_NB(b_str,2,arr_sz);
			rem--;
			//b_str<<=2;
			break;


		case ('C'):case ('c'):case ('1'):
			b_str[beg_arr_idx]<<=2;//L_shift_NB(b_str,2,arr_sz);
			++b_str[beg_arr_idx];//++b_str[arr_sz-1];
			rem--;
//			++(b_str<<=2);
			break;


		case 'G':case 'g':case '2':
			b_str[beg_arr_idx]<<=1;//L_shift_NB(b_str,1,arr_sz);
			++b_str[beg_arr_idx];//++b_str[arr_sz-1];
			b_str[beg_arr_idx]<<=1;//L_shift_NB(b_str,1,arr_sz);
			rem--;//(++(b_str<<=1))<<=1;
			break;

		case 'T':case 't':case '3':
			b_str[beg_arr_idx]<<=1;//L_shift_NB(b_str,1,arr_sz);
			++b_str[beg_arr_idx];
			b_str[beg_arr_idx]<<=1;//L_shift_NB(b_str,1,arr_sz);
			++b_str[beg_arr_idx];
			rem--;
			//++((++(b_str<<=1))<<=1);
			break;
		default:
			return b_str;
		}

	//	cout<<b_str<<endl;
	}
	return b_str;
}
void Init_HT(struct hashtable* ht,size_t ht_sz)
{
	ht->ht_sz=ht_sz;
	ht->store_pos=(struct bucket**)calloc(ht_sz,sizeof(struct bucket*));

	for(size_t i=0;i<ht_sz;++i)
	{
		ht->store_pos[i]=NULL;
	}
}

bool look_up_in_a_list(uint64_t seq,struct bucket *** ptr)
{
	bool found=0;
	//struct bucket ** bktptr;


	while((**ptr)!=NULL)
	{
		if((**ptr)->kmer_t.kmer==seq)
		{
		//	bktptr=ptr;
			break;
		}

		(*ptr)=&((**ptr)->nxt_bucket);
	}
	if((**ptr)==NULL)
	{
		found=0;
		//bktptr=ptr;
	}
	else
	{
		found=1;
	}
	return found;
}



// sparse k-mer graph construction function, it contains 2 rounds. In the first round we select the sparse k-mers,
//in the second, we build the links between the k-mers.
void Kmer_HT_Construction(uint64_t kmer ,struct hashtable *ht,int64_t *bucket_count,int K_size)
{

	size_t ht_sz=ht->ht_sz;
	bool flip=0,found=0;
	size_t hash_idx;
	bucket ** bktptr;
	uint64_t f_kmer,hv;

	f_kmer=get_rev_comp_seq(kmer,K_size);
	if(kmer>f_kmer)
	{
		uint64_t t=kmer;
		kmer=f_kmer;
		f_kmer=t;
		flip=1;
	}


	hv=MurmurHash64A(&kmer,sizeof(kmer),0);

	hash_idx=(size_t) (hv%ht_sz);

	bktptr= &(ht->store_pos[hash_idx]);

	found=look_up_in_a_list(kmer, &bktptr);



	if(found==0)
	{

		*(bktptr)=(struct bucket*)malloc(sizeof(struct bucket));
		memset(*bktptr,0,sizeof(struct bucket));
		(*(bktptr))->nxt_bucket=NULL;
		(*(bktptr))->kmer_t.kmer=kmer;
		(*bucket_count)++;

	}

}

bool Kmer_HT_Search(uint64_t kmer ,struct hashtable *ht,int K_size)
{

	size_t ht_sz=ht->ht_sz;
	bool flip=0,found=0;
	size_t hash_idx;
	bucket ** bktptr;
	uint64_t f_kmer,hv;

	f_kmer=get_rev_comp_seq(kmer,K_size);
	if(kmer>f_kmer)
	{
		uint64_t t=kmer;
		kmer=f_kmer;
		f_kmer=t;
		flip=1;
	}


	hv=MurmurHash64A(&kmer,sizeof(kmer),0);

	hash_idx=(size_t) (hv%ht_sz);

	bktptr= &(ht->store_pos[hash_idx]);

	found=look_up_in_a_list(kmer, &bktptr);

	return found;

}


int main(int argc, char* argv[])
{

	cout<<"Command:"<<endl;
	cout<<"ProgramFile  -k INPUT_FILE1 -o OUTPUT_PREFIX -s SINGLE_END_FILE_INPUT -s1 QUERY_PAIR1 -s2 QUERY_PAIR2 -fa CONVERT_FQ2FA -MinMatch MINIMUM_MATCHING_KMERS_PER_READ"<<endl;
	cout<<"add -fa 1 if you want to output fa instead of fq."<<endl;
//
	bool FA=0;
	string inputs;
	vector<string> in_filenames,search_filenames1,search_filenames2,search_filenames;
	string dirc,out_filename;
	int MinMatch = 1;
	for(int i=1;i<argc;++i)
	{

		if(strcmp(argv[i],"-k")==0)
		{
			i++;
			in_filenames.push_back(argv[i]);
			continue;
		}
		if(strcmp(argv[i],"-o")==0)
		{
			i++;
			out_filename=argv[i];

		}
		if(strcmp(argv[i],"-s1")==0)
		{
			i++;
			search_filenames1.push_back(argv[i]);

		}
		if(strcmp(argv[i],"-s2")==0)
		{
			i++;
			search_filenames2.push_back(argv[i]);

		}
		if(strcmp(argv[i],"-s")==0)
		{
			i++;
			search_filenames.push_back(argv[i]);

		}
		if(strcmp(argv[i],"-fa")==0)
		{
			i++;
			FA=atoi(argv[i]);

			continue;
		}
		if (strcmp(argv[i], "-MinMatch") == 0)
		{
			i++;
			MinMatch = atoi(argv[i]);

			continue;
		}


	}

	if (in_filenames.empty()||out_filename.empty()||(search_filenames1.empty()&&search_filenames.empty()))
	{
		cout<<"Too few inputs!"<<endl;
		return -1;
	}

	string suffix;


	bool paired_mode=0;
	if(search_filenames.empty())
	{
		paired_mode=1;
		cout<<search_filenames1.size()<<" input file(s)."<<endl;
		if (search_filenames1[0][search_filenames1[0].size() - 1] == 'a' || search_filenames1[0][search_filenames1[0].size() - 1] == 'A' || FA == 1)
		{
			suffix = ".fa";
		}
		//strncpy(res, s + index, n + 1);
		if (search_filenames1[0].substr(search_filenames1[0].size() -4, 4) == "a.gz" || search_filenames1[0].substr(search_filenames1[0].size() -4, 4) == "A.gz" || FA == 1)
		{
			suffix = ".fa.gz";
		}
		if ((search_filenames1[0][search_filenames1[0].size() - 1] == 'q' || search_filenames1[0][search_filenames1[0].size() - 1] == 'Q' )&& FA == 0)
		{
			suffix = ".fq";
		}
		if (search_filenames1[0].substr(search_filenames1[0].size() -4, 4) == "q.gz" || search_filenames1[0].substr(search_filenames1[0].size() -4, 4) == "Q.gz" || FA == 0)
		{
			suffix = ".fq.gz";
		}
	}
	else
	{
		paired_mode=0;
		cout<<search_filenames.size()<<" inputs."<<endl;
		if (search_filenames[0][search_filenames[0].size() - 1] == 'a' || search_filenames[0][search_filenames[0].size() - 1] == 'A' || FA == 1)
		{
			suffix = ".fa";
		}
		if (search_filenames[0].substr(search_filenames[0].size() -4, 4) == "a.gz" || search_filenames[0].substr(search_filenames[0].size() -4, 4) == "A.gz" || FA == 1)
		{
			suffix = ".fa.gz";
		}
		if (search_filenames[0][search_filenames[0].size() - 1] == 'q' || search_filenames[0][search_filenames[0].size() - 1] == 'Q' || FA == 0)
		{
			suffix = ".fq";
		}
		if (search_filenames[0].substr(search_filenames[0].size() -4, 4) == "q.gz" || search_filenames[0].substr(search_filenames[0].size() -4, 4) == "Q.gz" || FA == 0)
		{
			suffix = ".fq.gz";
		}
	}

	uint64_t *Kmer_bits;
	size_t sz_Kmer_bits=sizeof(uint64_t)*1000;
	Kmer_bits=(uint64_t *)malloc(sz_Kmer_bits);


	string in_line;
	char Kmer[Max_Kmer_Size];
	int len=0;
	uint64_t Kmer_cnt=0;
	int k=0;

	size_t curcnt=0;

	for (size_t in_cnt=0;in_cnt<in_filenames.size();++in_cnt)
	{
		//ifstream infile1(in_filenames[in_cnt].c_str());
		igzstream infile1(in_filenames[in_cnt].c_str());
		while(getline(infile1,in_line))
		{
			if(in_line[in_line.size()-1]=='\n'||in_line[in_line.size()-1]=='\r')
			{
				in_line.resize(in_line.size()-1);
			}
			size_t Read_sz=in_line.size();

			for (size_t i=0;i<Read_sz;++i)
			{
				while(i<Read_sz&&!isalpha(in_line[i]))
				{
					++i;
				}
				if (i==Read_sz)
				{break;}

				len=0;
				while(i<Read_sz&&isalpha(in_line[i]))
				{
					Kmer[len++]=in_line[i++];
				}

				Kmer[len]='\0';

				Kmer_cnt++;

			}
		}
		infile1.close();

	}
	int K_size=len;
	cout<<Kmer_cnt<<" Kmers."<<endl;
	cout<<"Kmer size: "<<K_size<<endl;

	struct hashtable ht;

	Init_HT(&ht,Kmer_cnt*6/5);

	uint64_t KmerBits;
	int64_t bucket_count=0;


	for (size_t in_cnt=0;in_cnt<in_filenames.size();++in_cnt)
	{
		//ifstream infile1(in_filenames[in_cnt].c_str());
		igzstream infile1(in_filenames[in_cnt].c_str());
		while(getline(infile1,in_line))
		{
			if(in_line[in_line.size()-1]=='\n'||in_line[in_line.size()-1]=='\r')
			{
				in_line.resize(in_line.size()-1);
			}
			size_t Read_sz=in_line.size();

			for (size_t i=0;i<Read_sz;++i)
			{
				while(i<Read_sz&&!isalpha(in_line[i]))
				{
					++i;
				}
				if (i==Read_sz)
				{break;}

				len=0;
				while(i<Read_sz&&isalpha(in_line[i]))
				{
					Kmer[len++]=in_line[i++];
				}

				Kmer[len]='\0';

				str2bitsarr(Kmer,K_size,&KmerBits,1);

				Kmer_HT_Construction(KmerBits,&ht,&bucket_count, K_size);

			}
		}
		infile1.close();

	}



	cout<<bucket_count<<" different Kmers."<<endl;

	if(paired_mode)
	{
		string in_line2,in_seq1,in_seq2;

		size_t Reads_found=0;
		string out_filename1,out_filename2;
		out_filename1 = out_filename + "_R1" + suffix;
		out_filename2 = out_filename + "_R2" + suffix;
		//ofstream outfile1(out_filename1.c_str());
		ogzstream outfile1(out_filename1.c_str());
		//ofstream outfile2(out_filename2.c_str());
		ogzstream outfile2(out_filename2.c_str());
		//int nr=0;
		for (size_t s_cnt=0;s_cnt<search_filenames1.size();++s_cnt)
		{
			cout<<"Processing file:" <<s_cnt+1<<endl;
			//ifstream infile2(search_filenames1[s_cnt].c_str());
			igzstream infile2(search_filenames1[s_cnt].c_str());
			//ifstream infile3(search_filenames2[s_cnt].c_str());
			igzstream infile3(search_filenames2[s_cnt].c_str());
			while(getline(infile2,in_line))
			{
				if(in_line[in_line.size()-1]=='\n'||in_line[in_line.size()-1]=='\r')
				{
					in_line.resize(in_line.size()-1);
				}
				getline(infile3,in_line2);

				if(in_line2[in_line2.size()-1]=='\n'||in_line2[in_line2.size()-1]=='\r')
				{
					in_line2.resize(in_line2.size()-1);
				}
				if (in_line[0]!='@'&&in_line[0]!='>')
				{
					continue;
				}
				else
				{
					getline(infile2,in_seq1);
					if(in_seq1[in_seq1.size()-1]=='\n'||in_seq1[in_seq1.size()-1]=='\r')
					{
						in_seq1.resize(in_seq1.size()-1);
					}
					getline(infile3,in_seq2);
					if(in_seq2[in_seq2.size()-1]=='\n'||in_seq2[in_seq2.size()-1]=='\r')
					{
						in_seq2.resize(in_seq2.size()-1);
					}
					uint64_t seq1,seq2;
					bool found=0,found1=0,found2=0;
				//	nr++;
					int seq1_sz=in_seq1.size();
					for(int k=0;k<seq1_sz+1-K_size;++k)
					{
						string KmerStr1=in_seq1.substr(k,K_size);
						bool Nflag=0;
						for(int kk=0;kk<KmerStr1.size();++kk)
						{
							if(KmerStr1[kk]=='N')
							{

								Nflag=1;
								break;
							}
						}
						if(Nflag==1)
						{continue;}

						str2bitsarr(KmerStr1.c_str(),K_size,&seq1,1);

						found1=Kmer_HT_Search(seq1 ,&ht,K_size);

						if (found1==1)
						{
							break;

						}

					}
					if (found1==0||MinMatch>1)
					{

						int seq2_sz=in_seq1.size();
						for(int k=0;k<seq2_sz-K_size+1;++k)
						{

							string KmerStr2=in_seq2.substr(k,K_size);

							bool Nflag=0;
							for(int kk=0;kk<KmerStr2.size();++kk)
							{
								if(KmerStr2[kk]=='N')
								{

									Nflag=1;
									break;
								}
							}
							if(Nflag==1)
							{continue;}



							str2bitsarr(KmerStr2.c_str(),K_size,&seq2,1);

							found2=Kmer_HT_Search(seq2 ,&ht,K_size);

							if (found2==1)
							{
								break;

							}
						}

					}
					if (MinMatch >= 2)
					{
						found = found1 & found2;

					}
					else
					{
						found = found1 | found2;
					}

					if (found==1)
					{
						//cout<<(int)found1<<" "<<(int)found2<<" "<<nr<<endl;
						if(FA)
						{
							in_line[0]='>';
						}
						outfile1<<in_line<<endl<<in_seq1<<endl;
						outfile2<<in_line<<endl<<in_seq2<<endl;
						if(in_line[0]=='@'&&(!FA))
						{
							getline(infile2,in_seq1);
							if(in_seq1[in_seq1.size()-1]=='\n'||in_seq1[in_seq1.size()-1]=='\r')
							{
								in_seq1.resize(in_seq1.size()-1);
							}
							outfile1<<in_seq1<<endl;
							getline(infile2,in_seq1);
							if(in_seq1[in_seq1.size()-1]=='\n'||in_seq1[in_seq1.size()-1]=='\r')
							{
								in_seq1.resize(in_seq1.size()-1);
							}
							outfile1<<in_seq1<<endl;

							getline(infile3,in_seq2);
							if(in_seq2[in_seq2.size()-1]=='\n'||in_seq2[in_seq2.size()-1]=='\r')
							{
								in_seq2.resize(in_seq2.size()-1);
							}
							outfile2<<in_seq2<<endl;
							getline(infile3,in_seq2);
							if(in_seq2[in_seq2.size()-1]=='\n'||in_seq2[in_seq2.size()-1]=='\r')
							{
								in_seq2.resize(in_seq2.size()-1);
							}
							outfile2<<in_seq2<<endl;
						}
						Reads_found++;

					}
				}

			//	size_t Read_sz=in_line.size();

			}
		}
		cout<<Reads_found<<" reads found."<<endl;
	}
	else
	{
		string in_line2,in_seq1;

		size_t Reads_found=0;
		string out_filename1;
		out_filename1 = out_filename + suffix;
		//ofstream outfile1(out_filename1.c_str());
		ogzstream outfile1(out_filename1.c_str());
		//int nr=0;

		for (size_t s_cnt=0;s_cnt<search_filenames.size();++s_cnt)
		{
			cout<<"Processing file:" <<s_cnt+1<<endl;
			//ifstream infile2(search_filenames[s_cnt].c_str());
			igzstream infile2(search_filenames[s_cnt].c_str());
			while(getline(infile2,in_line))
			{



				if (in_line[0]!='@'&&in_line[0]!='>')
				{
					continue;
				}
				else
				{

					getline(infile2,in_seq1);
					if(in_seq1[in_seq1.size()-1]=='\n'||in_seq1[in_seq1.size()-1]=='\r')
					{
						in_seq1.resize(in_seq1.size()-1);
					}

					uint64_t seq1,seq2;
					bool found=0,found1=0,found2=0;
			//		nr++;
					int n_kmer_matches = 0;
					int last_match = -10;
					int seq1_sz=in_seq1.size();
					for(int k=0;k<seq1_sz+1- K_size;++k)
					{
						string KmerStr1=in_seq1.substr(k,K_size);
						bool Nflag=0;
						for(int kk=0;kk<KmerStr1.size();++kk)
						{
							if(KmerStr1[kk]=='N')
							{

								Nflag=1;
								break;
							}
						}
						if(Nflag==1)
						{continue;}


						str2bitsarr(KmerStr1.c_str(),K_size,&seq1,1);

						found1=Kmer_HT_Search(seq1 ,&ht,K_size);

						if (found1==1)
						{
							if (last_match!=k-1)
							{
								n_kmer_matches++;
								last_match = k;
							}

						}



					}
					if (n_kmer_matches>=MinMatch)
					{
						found =1;
					}
					else
					{
						found = 0;
					}

					if (found==1)
					{
						if(FA)
						{
							in_line[0]='>';
						}

						//cout<<(int)found1<<" "<<" "<<nr<<endl;
						outfile1<<in_line<<endl<<in_seq1<<endl;



						if(in_line[0]=='@'&&(!FA))
						{
							getline(infile2,in_seq1);
							if(in_seq1[in_seq1.size()-1]=='\n'||in_seq1[in_seq1.size()-1]=='\r')
							{
								in_seq1.resize(in_seq1.size()-1);
							}
							outfile1<<in_seq1<<endl;
							getline(infile2,in_seq1);
							if(in_seq1[in_seq1.size()-1]=='\n'||in_seq1[in_seq1.size()-1]=='\r')
							{
								in_seq1.resize(in_seq1.size()-1);
							}
							outfile1<<in_seq1<<endl;
						}
						Reads_found++;

					}
				}

			//	size_t Read_sz=in_line.size();

			}
		}
		cout<<Reads_found<<" reads found."<<endl;
	}

	return 0;
}
