#include "residualjp2.hh"
#include "ycbcr.hh"
#include "ppm.hh"
#include "fileaux.hh"
#include "clip.hh"
#include "medianfilter.hh"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <vector>

#define CLEVELS 6
#define USE_JP2_DICTIONARY 0

void getJP2Header(unsigned char *JP2, unsigned char *&header, int JP2Size, int &headerSize) {

	for (int ii = 0; ii < JP2Size-1; ii++) {
		if ((unsigned short)(JP2[ii] << 8 | JP2[ii + 1]) == 0xFF90) { /*we have first tile*/
			headerSize = ii + 2;
			header = new unsigned char[headerSize];
			memcpy(header, JP2, headerSize);
			return;
		}
	}

	return;
}

int getJP2DictionaryIndex(unsigned char *header, int headerSize, 
	std::vector< std::vector<unsigned char>> JP2_dict) {

	for (int ii = 0; ii < JP2_dict.size(); ii++) {

		int  L = (int)JP2_dict.at(ii).size();

		if (L == headerSize) {
			if (memcmp(header, &JP2_dict.at(ii)[0], L) == 0) {
				return  ii;
			}
		}

	}

	return -1;
}

void updateJP2Dictionary(std::vector< std::vector<unsigned char>> &JP2_dict, unsigned char *JP2header, int headerSize) {

	std::vector< unsigned char > new_dict_element;

	for (int ii = 0; ii < headerSize; ii++) {
		new_dict_element.push_back(JP2header[ii]);
	}

	JP2_dict.push_back(new_dict_element);

}

void readResidualFromDisk(const char *jp2_residual_path_jp2, int &n_bytes_residual, FILE *input_LF, 
	std::vector< std::vector<unsigned char>> &JP2_dict) {

	int n_bytes_JP2 = 0;
	unsigned char *jp2_residual = 0;

	if (USE_JP2_DICTIONARY) {

		int dict_index = 0, headerSize = 0;
		unsigned char dict_index_char = 0;

		n_bytes_residual += (int)fread(&dict_index_char, sizeof(unsigned char), 1, input_LF) * sizeof(unsigned char);

		dict_index = (int)dict_index_char;

		if (JP2_dict.size() == 0 || dict_index > (int)JP2_dict.size()-1 ) {

			n_bytes_residual += (int)fread(&headerSize, sizeof(int), 1, input_LF) * sizeof(int);

			unsigned char *JP2header = new unsigned char[headerSize]();

			n_bytes_residual += (int)fread(JP2header, sizeof(unsigned char), headerSize, input_LF) * sizeof(unsigned char);

			updateJP2Dictionary(JP2_dict, JP2header, headerSize);

			delete[](JP2header);

		}

		unsigned char *jp2_residual_tmp;

		n_bytes_residual += (int)fread(&n_bytes_JP2, sizeof(int), 1, input_LF) * sizeof(int);
		jp2_residual_tmp = new unsigned char[n_bytes_JP2]();
		n_bytes_residual += (int)fread(jp2_residual_tmp, sizeof(unsigned char), n_bytes_JP2, input_LF);

		headerSize = (int)JP2_dict.at(dict_index).size();
		jp2_residual = new unsigned char[n_bytes_JP2 + headerSize]();

		memcpy(jp2_residual, &JP2_dict.at(dict_index)[0], headerSize);
		memcpy(jp2_residual + headerSize, jp2_residual_tmp, n_bytes_JP2);

		n_bytes_JP2 += headerSize;

		delete[](jp2_residual_tmp);

	}
	else {

		n_bytes_residual += (int)fread(&n_bytes_JP2, sizeof(int), 1, input_LF) * sizeof(int);
		jp2_residual = new unsigned char[n_bytes_JP2]();
		n_bytes_residual += (int)fread(jp2_residual, sizeof(unsigned char), n_bytes_JP2, input_LF)* sizeof(unsigned char);

	}

	FILE *jp2_res_file;
	jp2_res_file = fopen(jp2_residual_path_jp2, "wb");
	fwrite(jp2_residual, sizeof(unsigned char), n_bytes_JP2, jp2_res_file);
	fclose(jp2_res_file);

	delete[](jp2_residual);
}

void writeResidualToDisk(const char *jp2_residual_path_jp2, FILE *output_LF_file, int &n_bytes_residual, 
	std::vector< std::vector<unsigned char>> &JP2_dict) {

	int n_bytes_JP2 = aux_GetFileSize(jp2_residual_path_jp2);

	unsigned char *jp2_residual = new unsigned char[n_bytes_JP2]();
	FILE *jp2_color_residual_file = fopen(jp2_residual_path_jp2, "rb");
	fread(jp2_residual, sizeof(unsigned char), n_bytes_JP2, jp2_color_residual_file);
	fclose(jp2_color_residual_file);

	if (USE_JP2_DICTIONARY) {
		bool updateDictionary = false;
		/* get header */
		unsigned char *JP2header;
		int headerSize = 0;
		getJP2Header(jp2_residual, JP2header, n_bytes_JP2, headerSize);
		/* get index in dictionary */
		int dict_index = getJP2DictionaryIndex(JP2header, headerSize, JP2_dict);
		/* update dictionary if needed */
		if (dict_index == -1) {
			updateDictionary = true;
			updateJP2Dictionary(JP2_dict, JP2header, headerSize);
			dict_index = (int)JP2_dict.size() - 1;
			//printf("Dictonary update, index=%i\n", dict_index);
		}

		//printf("Using dictionary index:\t%i\n", dict_index);

		delete[](JP2header);

		/* write index of dictionary to bitstream */
		unsigned char dict_index_char = (unsigned char)dict_index;
		n_bytes_residual += (int)fwrite(&dict_index_char, sizeof(unsigned char), 1, output_LF_file) * sizeof(unsigned char);

		if (updateDictionary) { /* add update information if necessary */
			n_bytes_residual += (int)fwrite(&headerSize, sizeof(int), 1, output_LF_file) * sizeof(int);
			n_bytes_residual += (int)fwrite(&JP2_dict.at(dict_index)[0], sizeof(unsigned char), headerSize, output_LF_file) * sizeof(unsigned char);
		}

		/* write to codestream without header */

		n_bytes_JP2 = n_bytes_JP2 - headerSize;

		n_bytes_residual += (int)fwrite(&n_bytes_JP2, sizeof(int), 1, output_LF_file) * sizeof(int);
		n_bytes_residual += (int)fwrite(jp2_residual+ headerSize, sizeof(unsigned char), n_bytes_JP2, output_LF_file) * sizeof(unsigned char);
	
	}
	else {

		n_bytes_residual += (int)fwrite(&n_bytes_JP2, sizeof(int), 1, output_LF_file) * sizeof(int);
		n_bytes_residual += (int)fwrite(jp2_residual, sizeof(unsigned char), n_bytes_JP2, output_LF_file) * sizeof(unsigned char);


		/*n_bytes_JP2 = n_bytes_JP2 - jp2_headersize;

		n_bytes_residual += (int)fwrite(&n_bytes_JP2, sizeof(int), 1, output_LF_file) * sizeof(int);
		n_bytes_residual += (int)fwrite(jp2_residual + jp2_headersize, sizeof(unsigned char), n_bytes_JP2, output_LF_file) * sizeof(unsigned char);*/

	}

	delete[](jp2_residual);
	
}

void decodeResidualJP2(unsigned short *ps, const char *kdu_expand_path, const char *jp2_residual_path_jp2, const char *ppm_residual_path, int ncomp, const int offset, const int maxvali,
	const bool RESIDUAL_16BIT_bool)
{
	/* decode residual with kakadu */
	char kdu_expand_s[1024];
	sprintf(kdu_expand_s, "\"%s\"%s%s%s%s", kdu_expand_path, " -i ", jp2_residual_path_jp2, " -o ", ppm_residual_path);

	//std::cout << kdu_expand_s << "\n";

	int status = system_1(kdu_expand_s);

	signed int dv = RESIDUAL_16BIT_bool ? 1 : 2;
	signed int BP = RESIDUAL_16BIT_bool ? 16 : 10;
	signed int maxval = (1 << BP) - 1;// pow(2, BP) - 1;

	/* apply residual */

	unsigned short* jp2_residual;

	int nc1, nr1;

	if (aux_read16PGMPPM(ppm_residual_path, nc1, nr1, ncomp, jp2_residual))
	{

		for (int iir = 0; iir < nc1*nr1 * ncomp; iir++)
		{
			signed int val = (signed int)*(ps + iir) + (signed int)(jp2_residual[iir] * dv) - offset; // we assume that for 10bit case we have offset as 2^10-1, so go from 2^11 range to 2^10 and lose 1 bit of precision
			val = clip(val, 0, maxvali);
			*(ps + iir) = (unsigned short)(val);
		}

		delete[](jp2_residual);
	}
}

void decodeResidualJP2_YUV(unsigned short *ps, const char *kdu_expand_path, char *ycbcr_jp2_names[], char *ycbcr_pgm_names[], const int ncomp, const int offset, const int maxvali,
	const bool RESIDUAL_16BIT_bool)
{
	/* decode residual with kakadu */
	char kdu_expand_s[1024];

	for (int icomp = 0; icomp < ncomp; icomp++) {
		sprintf(kdu_expand_s, "\"%s\"%s%s%s%s", kdu_expand_path, " -i ", ycbcr_jp2_names[icomp], " -o ", ycbcr_pgm_names[icomp]);
		int status = system_1(kdu_expand_s);
		if (status < 0) {
			printf("KAKADU ERROR\nTERMINATING ... \n");
			exit(0);
		}
	}

	unsigned short *ycbcr = NULL;

	unsigned short *jp2_residual;

	int nc1, nr1,ncomp1;

	signed int BP = RESIDUAL_16BIT_bool ? 16 : 10;

	for (int icomp = 0; icomp < ncomp; icomp++) {
		if (aux_read16PGMPPM(ycbcr_pgm_names[icomp], nc1, nr1, ncomp1, jp2_residual))
		{
			if (ycbcr == NULL) {
				ycbcr = new unsigned short[nc1*nr1*3]();
				for (int ii = nc1*nr1; ii < 3 * nc1*nr1; ii++) {
					*(ycbcr + ii) = 128 * (1 << (BP - 8)); /* initialize chrominance */
				}
			}

			if (YUV_422) {
				if (icomp < 1) {
					memcpy(ycbcr + icomp*nr1*nc1, jp2_residual, sizeof(unsigned short)*nr1*nc1);
				}
				else {
					for (int cc = 0; cc < nc1; cc++) {
						memcpy(ycbcr + cc * 2 * nr1 + icomp*nr1*nc1*2, jp2_residual + cc*nr1, sizeof(unsigned short)*nr1);
						memcpy(ycbcr + (cc*2+1)*nr1 + icomp*nr1*nc1*2, jp2_residual + cc*nr1, sizeof(unsigned short)*nr1);
					}
					nc1 = nc1 * 2; // since we keep using this variable...
				}
			}
			else {
				memcpy(ycbcr + icomp*nr1*nc1, jp2_residual, sizeof(unsigned short)*nr1*nc1);
			}

			

			delete[](jp2_residual);
		}
	}

	signed int dv = RESIDUAL_16BIT_bool ? 1 : 2;
	
	signed int maxval = (1 << (BP)) - 1;// pow(2, BP) - 1;

	//for (int ii = 0; ii < nr1*nc1*ncomp; ii++) {
	//	*(ycbcr + ii) = clip(*(ycbcr + ii), (unsigned short)0, (unsigned short)maxval);
	//}

	unsigned short *rgb = new unsigned short[nr1*nc1*3]();


	if (RESIDUAL_16BIT_bool) {
		YCbCr2RGB(ycbcr, rgb, nr1, nc1, 16);
	}
	else {
		YCbCr2RGB(ycbcr, rgb, nr1, nc1, 10);
	}

	/* apply residual */

	for (int iir = 0; iir < nc1*nr1 * ncomp; iir++)
	{
		signed int val = (signed int)*(ps + iir) + (signed int)(rgb[iir]*dv) - offset; // we assume that for 10bit case we have offset as 2^10-1, so go from 2^11 range to 2^10 and lose 1 bit of precision
		val = clip(val, 0, maxvali);
		*(ps + iir) = (unsigned short)(val);
	}

	delete[](ycbcr);
	delete[](rgb);
}


void encodeResidualJP2_YUV(const int nr, const int nc, unsigned short *original_intermediate_view, unsigned short *ps, char *ycbcr_pgm_names[],
	const char *kdu_compress_path, char *ycbcr_jp2_names[], const float residual_rate, const int ncomp, const int offset, float rate_a,
	const bool RESIDUAL_16BIT_bool)
{
	/*establish residual*/
	unsigned short *residual_image = new unsigned short[nr*nc * 3]();

	signed int dv = RESIDUAL_16BIT_bool ? 1 : 2;
	signed int BP = RESIDUAL_16BIT_bool ? 16 : 10;
	signed int maxval = (1 << BP) - 1;// pow(2, BP) - 1;

	for (int iir = 0; iir < nr*nc*3; iir++) {
		signed int res_val = ( (((signed int)*(original_intermediate_view + iir)) - ((signed int)*(ps + iir)) + offset) )/dv;
		res_val = clip(res_val, 0, maxval);
		*(residual_image + iir) = (unsigned short)(res_val);
	}

	unsigned short *ycbcr = new unsigned short[nr*nc*3]();

	if (RESIDUAL_16BIT_bool) {
		RGB2YCbCr(residual_image, ycbcr, nr, nc, 16);
	}
	else {
		RGB2YCbCr(residual_image, ycbcr, nr, nc, 10);
	}

	unsigned short *tmp_im = new unsigned short[nr*nc]();

	for (int icomp = 0; icomp < ncomp; icomp++) {

		float rateR = residual_rate;

		if (icomp < 1) {
			rateR = rate_a*rateR;
		}
		else {
			rateR = (1-rate_a)*rateR/2;
		}

		if (YUV_422) {
			if (icomp < 1) {
				memcpy(tmp_im, ycbcr + icomp*nr*nc, sizeof(unsigned short)*nr*nc);
				aux_write16PGMPPM(ycbcr_pgm_names[icomp], nc, nr, 1, tmp_im);
			}
			else {
				for (int cc = 0; cc < nc; cc += 2) {
					memcpy(tmp_im + (cc / 2)*nr, ycbcr + cc*nr + nr*nc*icomp, sizeof(unsigned short)*nr);
				}
				aux_write16PGMPPM(ycbcr_pgm_names[icomp], nc/2, nr, 1, tmp_im);
			}
		}
		else {
			memcpy(tmp_im, ycbcr + icomp*nr*nc, sizeof(unsigned short)*nr*nc);
			aux_write16PGMPPM(ycbcr_pgm_names[icomp], nc, nr, 1, tmp_im);
		}

		char kdu_compress_s[1024];
	
		sprintf(kdu_compress_s, "\"%s\"%s%s%s%s%s%f%s%d", kdu_compress_path, " -i ", ycbcr_pgm_names[icomp], " -o ", ycbcr_jp2_names[icomp], " -num_threads 0 -no_weights -full -no_info -precise -rate ", rateR,
			" Clevels=",CLEVELS);

		int status = system_1(kdu_compress_s);

		if (status < 0) {
			printf("KAKADU ERROR\nTERMINATING ... \n");
			exit(0);
		}


	}

	delete[](tmp_im);
	delete[](ycbcr);
	delete[](residual_image);


}

void encodeResidualJP2(const int nr, const int nc, unsigned short *original_intermediate_view, unsigned short *ps, const char *ppm_residual_path,
	const char *kdu_compress_path, const char *jp2_residual_path_jp2, const float residual_rate, const int ncomp, const int offset, const bool RESIDUAL_16BIT_bool)
{
	/*establish residual*/
	unsigned short *residual_image = new unsigned short[nr*nc * ncomp]();

	signed int dv = RESIDUAL_16BIT_bool ? 1 : 2;
	signed int BP = RESIDUAL_16BIT_bool ? 16 : 10;
	signed int maxval = (1 << BP) - 1;// pow(2, BP) - 1;

	for (int iir = 0; iir < nr*nc*ncomp; iir++) {
		signed int res_val = ((((signed int)*(original_intermediate_view + iir)) - ((signed int)*(ps + iir)) + offset)) / dv;
		res_val = clip(res_val, 0, maxval);
		*(residual_image + iir) = (unsigned short)(res_val);
	}

	aux_write16PGMPPM(ppm_residual_path, nc, nr, ncomp, residual_image);

	delete[](residual_image);

	/* here encode residual with kakadu */

	char kdu_compress_s[1024]; // tolerance 0 ?
	sprintf(kdu_compress_s, "\"%s\"%s%s%s%s%s%f%s%d", kdu_compress_path, " -i ", ppm_residual_path, " -o ", jp2_residual_path_jp2, " -no_weights -full -no_info -precise -rate ", residual_rate,
		" Clevels=", CLEVELS);

	//std::cout << kdu_compress_s << "\n";

	int status = system_1(kdu_compress_s);
}


void encodeResidualHM(
    const int nr, 
    const int nc, 
    unsigned short *original_intermediate_view, 
    unsigned short *ps, 
    const char *ppm_residual_path,
    const char *kdu_compress_path, 
    const char *jp2_residual_path_jp2, 
    const float residual_rate, 
    const int ncomp, 
    const int offset, 
    const bool RESIDUAL_16BIT_bool,
    const int Q)
{
    /*establish residual*/
    unsigned short *residual_image = new unsigned short[nr*nc * ncomp]();

    signed int dv = RESIDUAL_16BIT_bool ? 1 : 2;
    signed int BP = RESIDUAL_16BIT_bool ? 16 : 10;
    signed int maxval = (1 << BP) - 1;// pow(2, BP) - 1;

    dv = static_cast<signed int>(Q);

    for (int iir = 0; iir < nr*nc*ncomp; iir++) {
        signed int res_val = ((((signed int)*(original_intermediate_view + iir)) - ((signed int)*(ps + iir)) + offset)) / dv;
        res_val = clip(res_val, 0, maxval);
        *(residual_image + iir) = (unsigned short)(res_val);
    }

    int CUsize = 64;
    int nc1 = nc%CUsize>0 ? (nc / CUsize+ 1)*CUsize : (nc / CUsize )*CUsize;
    int nr1 = nr%CUsize>0 ? (nr / CUsize + 1)*CUsize : (nr / CUsize)*CUsize;

    unsigned short *temp_im = 
        new unsigned short[nc1*nr1*ncomp]();

    for (int ii = 0; ii < nr1*nc1*ncomp; ii++) {
        temp_im[ii] = offset / dv;
    }

    for (int rr = 0; rr < nr; rr++) {
        for (int cc = 0; cc < nc; cc++) {
            for (int icomp = 0; icomp < ncomp; icomp++) {
                temp_im[rr + (nr1)*cc + icomp*(nc1)*(nr1)] =
                    residual_image[rr + nr*cc + icomp*nc*nr];
            }
        }
    }

    //aux_write16PGMPPM(ppm_residual_path, nc, nr, ncomp, residual_image);
    aux_write16PGMPPM(ppm_residual_path, nc1, nr1, ncomp, temp_im);

    delete[](temp_im);
    delete[](residual_image);

    /* convert to YUV420 */

    char RGB_2_YUV420_s[1024];
    sprintf(
        RGB_2_YUV420_s,
        "%s%s%s",
        "C:/Local/astolap/Programs/ffmpeg-20190225-f948082-win64-static/bin/ffmpeg.exe -y -i ",
        ppm_residual_path,
        " -c:v rawvideo -pix_fmt yuv420p10le C:/Temp/tmp.yuv");

    int status = system_1(RGB_2_YUV420_s);

    /* make HM cfg */

    std::ifstream inFile("C:/Local/astolap/Data/JPEG_PLENO_2019/GENEVA/Matlab/intra_cfg_1.txt");
    std::ostringstream buffer;
    buffer << inFile.rdbuf();
    inFile.close();

    std::ofstream outFile("C:/Temp/intra_cfg_enc.cfg", std::ios_base::out);
    std::string str;
    str.append(buffer.str());
    str.append("\nTargetBitrate\t: ");
    str.append(std::to_string(static_cast<uint32_t>(residual_rate*nr*nc)));
    outFile << str;
    //outFile << "TargetBitrate\t: ";
    //outFile << static_cast<uint32_t>(residual_rate*nr*nc);
    outFile.close();

    /* encode with HM */

    char HM_call_s[1024];
    sprintf(HM_call_s,
        "C:/Local/astolap/Data/JPEG_PLENO_2019/BRUSSELS/HEVC-HM/bin/vc2015/x64/Release/TAppEncoder.exe -c C:/Temp/intra_cfg_enc.cfg -i C:/Temp/tmp.yuv -o C:/Temp/tmp_rec.yuv -fr 1 -wdt %i -hgt %i -b %s",
        nc1,
        nr1,
        jp2_residual_path_jp2);

    status = system_1(HM_call_s);

}

void decodeResidualHM(
    const int nr,
    const int nc,
    unsigned short *ps, 
    const char *kdu_expand_path, 
    const char *jp2_residual_path_jp2,
    const char *ppm_residual_path,
    int ncomp, 
    const int offset,
    const int maxvali,
    const bool RESIDUAL_16BIT_bool,
    int Q)
{

    /*decode HM*/
    char hm_decode_s[1024];
    sprintf(hm_decode_s,
        "C:/Local/astolap/Data/JPEG_PLENO_2019/BRUSSELS/HEVC-HM/bin/vc2015/x64/Release/TAppDecoder.exe -b %s -o %s",
        jp2_residual_path_jp2,
        "C:/Temp/tmp_rec.yuv");

    int status = system_1(hm_decode_s);

    int CUsize = 64;
    int nc1 = nc%CUsize>0 ? (nc / CUsize + 1)*CUsize : (nc / CUsize)*CUsize;
    int nr1 = nr%CUsize>0 ? (nr / CUsize + 1)*CUsize : (nr / CUsize)*CUsize;

    /*convert to rgb*/
    char YUV420_2_RGB_s[1024];
    sprintf(
        YUV420_2_RGB_s,
        "C:/Local/astolap/Programs/ffmpeg-20190225-f948082-win64-static/bin/ffmpeg.exe -y -s:v %ix%i -pix_fmt yuv420p10le -r 1 -i %s -y %s",
        nc1,
        nr1,
        "C:/Temp/tmp_rec.yuv",
        ppm_residual_path);

    status = system_1(YUV420_2_RGB_s);

    /* convert */

    signed int dv = RESIDUAL_16BIT_bool ? 1 : 2;
    signed int BP = RESIDUAL_16BIT_bool ? 16 : 10;
    signed int maxval = (1 << BP) - 1;// pow(2, BP) - 1;

                                      /* apply residual */

    dv = static_cast<signed int>(Q);


    unsigned short *jp2_residual_padded, *jp2_residual;

    aux_read16PGMPPM(ppm_residual_path, nc1, nr1, ncomp, jp2_residual_padded);

    jp2_residual = new unsigned short[nr*nc*ncomp]();

    for (int rr = 0; rr < nr; rr++) {
        for (int cc = 0; cc < nc; cc++) {
            for (int icomp = 0; icomp < ncomp; icomp++) {
                jp2_residual[rr+cc*nr+nr*nc*icomp] = 
                    jp2_residual_padded[rr + cc*nr1 + nr1*nc1*icomp] >> 6;
            }
        }
    }

    delete[](jp2_residual_padded);

    aux_write16PGMPPM(ppm_residual_path, nc, nr, ncomp, jp2_residual);

    if (aux_read16PGMPPM(ppm_residual_path, nc1, nr1, ncomp, jp2_residual))
    {

        for (int iir = 0; iir < nc1*nr1 * ncomp; iir++)
        {
            signed int val = (signed int)*(ps + iir) + (signed int)(jp2_residual[iir] * dv) - offset; // we assume that for 10bit case we have offset as 2^10-1, so go from 2^11 range to 2^10 and lose 1 bit of precision
            val = clip(val, 0, maxvali);
            *(ps + iir) = (unsigned short)(val);
        }

        
    }
    delete[](jp2_residual);
}


void encodeMonochromeResidualHM(
    const int nr,
    const int nc,
    unsigned short *original_intermediate_view,
    unsigned short *ps,
    const char *ppm_residual_path,
    const char *kdu_compress_path,
    const char *jp2_residual_path_jp2,
    const float residual_rate,
    const int ncomp,
    const int offset,
    const bool RESIDUAL_16BIT_bool,
    const int Q)
{
    /*establish residual*/
    unsigned short *residual_image = new unsigned short[nr*nc * ncomp]();

    signed int dv = RESIDUAL_16BIT_bool ? 1 : 2;
    signed int BP = RESIDUAL_16BIT_bool ? 16 : 10;
    signed int maxval = (1 << BP) - 1;// pow(2, BP) - 1;

    dv = static_cast<signed int>(Q);

    for (int iir = 0; iir < nr*nc*ncomp; iir++) {
        signed int res_val = ((((signed int)*(original_intermediate_view + iir)) - ((signed int)*(ps + iir)) + offset)) / dv;
        res_val = clip(res_val, 0, maxval);
        *(residual_image + iir) = (unsigned short)(res_val);
    }

    int CUsize = 64;
    int nc1 = nc%CUsize>0 ? (nc / CUsize + 1)*CUsize : (nc / CUsize)*CUsize;
    int nr1 = nr%CUsize>0 ? (nr / CUsize + 1)*CUsize : (nr / CUsize)*CUsize;

    unsigned short *temp_im =
        new unsigned short[nc1*nr1*ncomp]();

    for (int ii = 0; ii < nr1*nc1*ncomp; ii++) {
        temp_im[ii] = offset / dv;
    }

    for (int rr = 0; rr < nr; rr++) {
        for (int cc = 0; cc < nc; cc++) {
            for (int icomp = 0; icomp < ncomp; icomp++) {
                temp_im[rr + (nr1)*cc + icomp*(nc1)*(nr1)] =
                    residual_image[rr + nr*cc + icomp*nc*nr];
            }
        }
    }

    FILE *temp_yuv = fopen("C:/Temp/tmp.yuv", "wb");
    fwrite(residual_image, sizeof(unsigned short), nr1*nc1,temp_yuv);
    fclose(temp_yuv);

    delete[](residual_image);
    delete[](temp_im);

    /* make HM cfg */

    std::ifstream inFile("C:/Local/astolap/Data/JPEG_PLENO_2019/GENEVA/Matlab/intra_cfg_400.cfg");
    std::ostringstream buffer;
    buffer << inFile.rdbuf();
    inFile.close();

    std::ofstream outFile("C:/Temp/intra_cfg_enc.cfg", std::ios_base::out);
    std::string str;
    str.append(buffer.str());
    str.append("\nTargetBitrate\t: ");
    str.append(std::to_string(static_cast<uint32_t>(residual_rate*nr*nc)));
    outFile << str;
    //outFile << "TargetBitrate\t: ";
    //outFile << static_cast<uint32_t>(residual_rate*nr*nc);
    outFile.close();

    /* encode with HM */

    char HM_call_s[1024];
    sprintf(HM_call_s,
        "C:/Local/astolap/Data/JPEG_PLENO_2019/BRUSSELS/HEVC-HM/bin/vc2015/x64/Release/TAppEncoder.exe -c C:/Temp/intra_cfg_enc.cfg -i C:/Temp/tmp.yuv -o C:/Temp/tmp_rec.yuv -fr 1 -wdt %i -hgt %i -b %s",
        nr1,
        nc1,
        jp2_residual_path_jp2);

    int status = system_1(HM_call_s);

}


void decodeMonochromeResidualHM(
    const int nr,
    const int nc,
    unsigned short *ps,
    const char *kdu_expand_path,
    const char *jp2_residual_path_jp2,
    const char *ppm_residual_path,
    const int ncomp,
    const int offset,
    const int maxvali,
    const bool RESIDUAL_16BIT_bool,
    int Q)
{

    /*decode HM*/
    char hm_decode_s[1024];
    sprintf(hm_decode_s,
        "C:/Local/astolap/Data/JPEG_PLENO_2019/BRUSSELS/HEVC-HM/bin/vc2015/x64/Release/TAppDecoder.exe -b %s -o %s",
        jp2_residual_path_jp2,
        "C:/Temp/tmp_rec.yuv");

    int status = system_1(hm_decode_s);

    int CUsize = 64;
    int nc1 = nc%CUsize>0 ? (nc / CUsize + 1)*CUsize : (nc / CUsize)*CUsize;
    int nr1 = nr%CUsize>0 ? (nr / CUsize + 1)*CUsize : (nr / CUsize)*CUsize;

    unsigned short *tmp_rec_im = new unsigned short[nr1*nc1]();

    FILE *tmp_yuv_rec = fopen("C:/Temp/tmp_rec.yuv", "rb");
    fread(tmp_rec_im, sizeof(unsigned short), nr1*nc1, tmp_yuv_rec);
    fclose(tmp_yuv_rec);


    /* convert */

    signed int dv = RESIDUAL_16BIT_bool ? 1 : 2;
    signed int BP = RESIDUAL_16BIT_bool ? 16 : 10;
    signed int maxval = (1 << BP) - 1;// pow(2, BP) - 1;

                                      /* apply residual */

    dv = static_cast<signed int>(Q);


    unsigned short *jp2_residual;

    //aux_read16PGMPPM(ppm_residual_path, nc1, nr1, ncomp, jp2_residual_padded);

    jp2_residual = new unsigned short[nr*nc*ncomp]();

    for (int rr = 0; rr < nr; rr++) {
        for (int cc = 0; cc < nc; cc++) {
            for (int icomp = 0; icomp < ncomp; icomp++) {
                jp2_residual[rr + cc*nr + nr*nc*icomp] =
                    tmp_rec_im[rr + cc*nr1 + nr1*nc1*icomp];
            }
        }
    }

    delete[](tmp_rec_im);

    aux_write16PGMPPM(ppm_residual_path, nc, nr, ncomp, jp2_residual);

    int ncomp1;

    if (aux_read16PGMPPM(ppm_residual_path, nc1, nr1, ncomp1, jp2_residual))
    {

        for (int iir = 0; iir < nc1*nr1 * ncomp; iir++)
        {
            signed int val = (signed int)*(ps + iir) + (signed int)(jp2_residual[iir] * dv) - offset; // we assume that for 10bit case we have offset as 2^10-1, so go from 2^11 range to 2^10 and lose 1 bit of precision
            val = clip(val, 0, maxvali);
            *(ps + iir) = (unsigned short)(val);
        }

       
    }

    delete[](jp2_residual);
}