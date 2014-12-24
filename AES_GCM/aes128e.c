#include <stdint.h>
#include "aes128e.h"

/* Multiplication by two in GF(2^8). Multiplication by three is xtime(a) ^ a */
#define xtime(a) ( ((a) & 0x80) ? (((a) << 1) ^ 0x1b) : ((a) << 1) )

/* The S-box table */
static const unsigned char sbox[256] = {
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5,
    0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0,
    0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc,
    0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a,
    0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0,
    0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b,
    0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85,
    0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5,
    0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17,
    0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88,
    0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c,
    0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9,
    0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6,
    0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e,
    0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94,
    0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68,
    0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 };

/* The round constant table (needed in KeyExpansion) */
static const unsigned char rcon[10] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 
    0x20, 0x40, 0x80, 0x1b, 0x36 };
// Performs a matrix multiplication of the column of the state with the given matrix in the galois field
void mixColumn(unsigned char *column)
{
	unsigned char temp[4];
	unsigned char i;
	for(i=0;i<4;i++) 
	{
		temp[i] = column[i];
	}
	column[0] = xtime(temp[0]) ^ temp[3] ^ temp[2] ^ xtime(temp[1])^temp[1];
	column[1] = xtime(temp[1]) ^ temp[0] ^ temp[3] ^ xtime(temp[2])^temp[2];
	column[2] = xtime(temp[2]) ^ temp[1] ^ temp[0] ^ xtime(temp[3])^temp[3];
	column[3] = xtime(temp[3]) ^ temp[2] ^ temp[1] ^ xtime(temp[0])^temp[0];
}

// Perform the matrix multiplication for each column in the state
void mixColumns(unsigned char *state) 
{ 
	int i, j; 
	unsigned char column[4]; 
	// Iterate over columns 
	for (i = 0; i < 4; i++) 
	{ 
		// Create one column
		for (j = 0; j < 4; j++) 
		{ 
			column[j] = state[(i*4)+j]; 
		} 
		//apply the mixColumn on one column 
		mixColumn(column); 
		//put the values back into the state
		for (j = 0; j < 4; j++) 
		{ 
			state[(i*4)+j] = column[j]; 
		} 
	} 
} 
// getSBoxValue method will return the sbox value associated with the hexadecimal number num
unsigned char getSBoxValue(unsigned char num) 
{ 
	return sbox[num]; 
} 

// In this method, the rows are rotated by the number of times specified.
void shiftRow(unsigned char *row, int num)
{
	int i,j;
	unsigned char temp;
	for(i=0;i<num;i++)
	{
		temp=row[0];
		//printf("temp %x ",temp);
		for(j=0;j<3;j++)
			row[j]=row[j+1];
		row[3]=temp;
	}	
}
//In this method, the second row is left rotated once,
//third row is left rotated twice,
// fourth row is left rotated thrice
// the first row is unchanged
void shiftRows(unsigned char *state)
{
	int i,j;
	unsigned char mat[4][4];
	//First create a matrix of 4 rows and 4 columns from the state which is stored in 1 D
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			mat[j][i] = state[i*4+j];
	// Shift the rows based on the row number
	for(i=0;i<4;i++)
	{
		shiftRow(mat[i],i);
	}
	//The matrix is converted back to the 1 D array and stored in state
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			state[i*4+j] = mat[j][i];
}
// AddRoundKey method, the key in that round is XORed with the state
void addRoundKey(unsigned char *state,unsigned char *key)
{
	int i;
	for(i=0;i<16;i++)
	{
		state[i]=state[i]^key[i];
	}
}
// This method performs the Key expansion
void key_expansion(unsigned char *cipher_key) 
{
	unsigned char temp[4];
	// length is 16 because the first cipher key provided
	unsigned char length = 16;
	unsigned char i,k=0;
	// We need 11 sets of sixteen bytes each for AES-128 
	while(length < 176) 
	{
		// Copy the temporary variable over from the last 4-byte block
		for(i = 0; i < 4; i++) 
		{
				temp[i] = cipher_key[i + length - 4];
				//printf("temp[i] %x ",temp[i]);
		}
		// Every four blocks (of four bytes) perform the operation		  
		if(length % 16 == 0) 
		{
		// Rotate the vector and get the sboxvalue of the associated values.
			shiftRow(temp,1);
			for(i = 0; i < 4; i++) 
				temp[i] = getSBoxValue(temp[i]);
			// The first value is XORed with the rcon value
			temp[0]=temp[0]^rcon[k];
			for(i = 0; i < 4; i++) 
			{
				//printf("temp[i] %x+%x ",temp[i],rcon[k]);
			}
			k++;
		}
		//Finally add the calculated value back to the cipher_key
		for(i = 0; i < 4; i++) 
		{
			cipher_key[length] = cipher_key[length - 16] ^ temp[i];
			//printf("cy[i] %x+%x ",cipher_key[length], temp[i]);
			length++;
		}
	}
}
//Performs the first operation in AES-128, where every value in the state is replaced by its associated sBoxValue
void subsitute_Bytes(unsigned char *state) 
{
	int i;
	for(i=0;i<16;i++)
	{
		state[i]=getSBoxValue(state[i]);
	}
}
//The cipher key is a 176 byte length array. We need 16 bytes for each round. 
// getRoundKey function returns a 16 byte key for each round
void getRoundKey(unsigned char *expanded_key, unsigned char *round_key) 
{ 
	int i; 
	for (i = 0; i < 16; i++) 
	{ 
		round_key[i] = expanded_key[i];
	} 
} 

/* Under the 16-byte key at k, encrypt the 16-byte plaintext at p and store it at c. */
void aes128e(unsigned char *c, const unsigned char *p, const unsigned char *k) {
	unsigned char expanded_key[176];
	unsigned char round_key[16]; 
	unsigned char state[16];
	int i,j;
	for(i=0;i<16;i++)
	{
		expanded_key[i] = k[i];
		state[i] = p[i];
	}
	//First generate the keys required for all 10 rounds
	key_expansion(expanded_key);
	//Get the key for the intial round
	getRoundKey(expanded_key,round_key);
	//Perform the add round key first
	addRoundKey(state,round_key);
	// Perform four functions 9 times.
	for(j=1;j<10;j++)
	{
		getRoundKey(expanded_key + 16*j,round_key); //  Gets the key required for that particular round
		subsitute_Bytes(state);
		shiftRows(state);
		mixColumns(state); 
		addRoundKey(state, round_key); 
	}
	// For the last round mix columns is not performed
	getRoundKey(expanded_key + 16*10,round_key);
	subsitute_Bytes(state); 
	shiftRows(state);
	addRoundKey(state, round_key); 
	// Finally the values are stored in array c (pointer)
	for(i=0;i<16;i++)
	{
		c[i] = state[i];
	}
}


