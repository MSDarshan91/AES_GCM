
#include "aes128gcm.h"
//This function is used perform the multiplication in galois field.
// Inputs are X and Y the output is stored in out.
void galois_multiply(unsigned char* X, unsigned char* Y, unsigned char *out)
{
  unsigned char V[16];
  unsigned char Z[16];
  //Initialization of Z and V
  for(int i=0;i<16;i++) 
  {
    Z[i]=0x00;
    V[i]=Y[i];
  }

  for (int i = 0; i < 16; i++) 
  {
		for (int j = 0; j < 8; j++) 
		{
			// If the first bit of input X is 1
			if (X[i] & 1<<(7- j)) 
			{
				for(int k=0;k<16;k++)
					 Z[k]=Z[k]^V[k];
			} 
			//If the last bit of V is one
			if (V[15] & 0x01) 
			{
				//Right shifting all the values of the array V
				for(int k=0;k<16;k++)
				{ 
					if (k!=0)
						if (V[15-k] & 0x01)
							V[16-k] = V[16-k] | 0x80;
					V[15-k]=V[15-k]>>1;
				 }
				 //XOR with R
			  	 V[0]=V[0]^0xe1;
			} 
			else 
			{
				//Right shifting all the values of the array V
				for(int k=0;k<16;k++)
				{ 
					if (k!=0)
						if (V[15-k] & 0x01)
							V[16-k] = V[16-k] | 0x80;
					V[15-k]=V[15-k]>>1;
			       }	
			}
		}
	}
    for(int i=0;i<16;i++)
	{
	//	Assigning all the final values of Z to out.
	  out[i] = Z[i]; 
	 // printf("%x ",out[i]);
	}
}
//Function to perform counter. X is the input, K is the key, ICB is the Initial Counter Block, len_X is the length of array X in bytes. Y is the output
void g_counter(unsigned char *X,const unsigned char *k,unsigned char *ICB, long len_X,unsigned char *Y)
{
	unsigned char CB[16];
	unsigned char encrypt_CB[16];
	int inc = 0;
	unsigned int i32count = 1;
	//If len_X is 0 then return empty string
	if(len_X > 0)
	{
		//Initialization of the Counter block
		for(int i=0;i<16;i++)
			CB[i] = ICB[i];
		for(int i=0;i<len_X;i++)
		{
			//Incrementing the integer value to be appended to CB
			i32count++;
			if(i32count==0)
				i32count=1;
            CB[12]=i32count>>24;
            CB[13]=i32count>>16;
            CB[14]=i32count>>8;
            CB[15]=i32count;
			//Encrypt the CB
			aes128e(encrypt_CB,CB,k);
			for(int j =0;j<16;j++)
			{
				//Xor input with encrypted CB
				Y[(i*16)+j] = encrypt_CB[j] ^ X[(i*16)+j];
			}
		}
	}
}
//This function is used to perform the hashing
// H is the Hash Subkey, X is the input.
// m is the length of X and Y_m is the output
void g_hash(unsigned char *H, unsigned char *X, unsigned int m,unsigned char *Y_m)
{
  unsigned char Y[16];
  //Initialization of Y
  for( int i=0;i<16;i++)
  {
		Y[i] = 0x00;
  }
  for(int i=0; i < m; i++)
  {
	unsigned char temp[16];
	//XOR input with Y
	for(int j =0;j<16;j++)
	{
		temp[j] = Y[j] ^ X[(i*16)+j];
	}
	//Perform the galois multiply on H and temp and store in Y.
    galois_multiply(temp, H, Y);
  }
  for( int i=0;i<16;i++)
  {
  //Finally assign the value of Y to output.
		Y_m[i] = Y[i];
		//printf("%x",Y_m[i]);
  }
//  printf("\n");
}

/* Under the 16-byte (128-bit) key "k", 
and the 12-byte (96-bit) initial value "IV", 
encrypt the plaintext "plaintext" and store it at "ciphertext". 
The length of the plaintext is a multiple of 16-byte (128-bit) given by len_p (e.g., len_p = 2 for a 32-byte plaintext). 
The length of the ciphertext "ciphertext" is len_p*16 bytes. 
The authentication tag is obtained by the 16-byte tag "tag". 
For the authentication an additional data "add_data" can be added. 
The number of blocks for this additional data is "len_ad" (e.g., len_ad = 1 for a 16-byte additional data). 
*/
void aes128gcm(unsigned char *ciphertext, unsigned char *tag, const unsigned char *k, const unsigned char *IV, const unsigned char *plaintext, const unsigned long len_p, const unsigned char* add_data, const unsigned long len_ad) 
{
 
	unsigned char J_0[16]; 
	unsigned char X[16*(len_p+len_ad+1)];
	//Initial Counter block J_0
	for(int i = 0;i<16;i++)
	{
		J_0[i] = 0x00;
		if(i<12)
			J_0[i] = IV[i];
	}
	J_0[15]=0x01;
	unsigned char H[16];
	unsigned char zero[16];
	for(int i = 0;i<16;i++)
	{
		zero[i]=0x00;
	}
	unsigned char P[16*len_p];
	for(int i=0; i<16*len_p; i++)
		P[i]= plaintext[i];
	//Encrypt zero vector with key k and store in H
	aes128e(H,zero, k);
	
	g_counter(P,k,J_0,len_p,ciphertext);
	//Appending cipher text and additional data to X
	for(int i=0; i<16*len_ad; i++)
		X[i]= add_data[i];
	for(int i=16*len_ad; i<16*(len_ad+len_p); i++)
		X[i]= ciphertext[i-(16*len_ad)];
	//Finally append the length of additional data and length of input data to X
	for(int i=16*(len_ad+len_p);i< (16*(len_ad+len_p))+7;i++)
	{
		X[i] = (16*8*len_ad) >> (8*(7-(i-16*(len_ad+len_p)))) & 0xff ;
	}
	X[(16*(len_ad+len_p))+7] = (16*8*len_ad)  & 0xff;
	for(int i=(16*(len_ad+len_p))+8;i< (16*(len_ad+len_p))+15;i++)
	{
		X[i] = (16*8*len_p) >> (8*(7-(i-(16*(len_ad+len_p)))+8)) & 0xff ;
	}
	X[(16*(len_ad+len_p))+15] = (16*8*len_p)  & 0xff;
	unsigned long len_x =(len_ad+len_p+1);
	unsigned char Y[16*3];
	for(int i=0;i<len_x*16;i++)
		//printf("%x",X[i]);
	//Perform the hashing on X
	g_hash(H, X,len_x, Y);
	//Performing the first stage of G Counter instead of calling the G counter function and getting the MSB of the output.
	unsigned char temp[16*3];
	unsigned char en_counter[16];
	aes128e(en_counter,J_0,k);
	for(int i=0;i<16;i++)
	{
		tag[i] = en_counter[i]^Y[i];
		//printf("%x",Y[i]);
	}
}



