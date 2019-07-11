//////////////////////////////////
////// Serebryakov S., 2019 //////
//////////////////////////////////
// compilation: gcc -lncurses -lm -lcfitsio focus_controller.c tat_info.c fli_func.c ccd_func.c -o focus_controller
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>
#include "fli_func.h"
#include "ccd_func.h"
#include "tat_info.h"
#include "symblc_const.h"
/////////////////////////////
#define ccd_limit_val   65535 
#define focuser_speed   0.0015 // seconds per step
////////////////////////
st_tat_info *p_tat_info;       
///////////////
typedef struct{
	float x_cr;
	float y_cr;
	float flux;
	int countp;
}blob_el; 
typedef struct{
	int    hist_len;
	int   *hist_cnt;
	float *hist_val;
	float *hist_cdf;
}histo;
////////////////////////////////////////			
float otsuThreshold(float*,float*,long);			 
float catchmentArea(float*,long*,blob_el*);
int blobDetector(float*,blob_el*,long*,int,int,int);	
int findallBlobs(float *,long*,blob_el**);
int extract_EAA(float*,long*,int,int,float);		
int subtrct_MBG(float*,long);
int equal_histo(float*,long);
float compute_HFD(float*,long*);
histo *histogram_maker(float*,long);
///////////////////////////////////////
int comp(const void *a, const void *b){
     return *(int*)a - *(int*)b;
}
///////////////////////////////// 
int main(int argc, char *argv[]){
	float exposure;
	long steps_per_request;
	if (argc != 4){
		printf("Example:%s def_expose steps_per_loop base_direction , where:\n",argv[0]);
		return EXIT_FAILURE;	
	} else {
		exposure = atof(argv[1]);	  
		steps_per_request = atof(argv[2]);   //500
		steps_per_request*=(atoi(argv[3]));
	}
	fitsfile *f_in,*f_out;
	long move_to, naxes[2] = {1,1}, fpx[2] = {1,1};
	int status_in=0,status_out=0, shmid, bitpix, naxis, binning=1,radius=100;
	int focused=0, saturated=0, index, ea_area, x_cr, y_cr, target_row, blob_count;
	float *img_arr, *curr_ptr, *bn_mask, *target_EAA, aperture = 1.0;
	float threshold, HFD[2] = {10,10}, arr_min =0,arr_max =0;
	char fli_cmd[20],sys_com[50],img_dir[40],img_nme[60];
	create_tat_info_shm( &shmid, &p_tat_info);
	//////////////////////////////////////////////
	if (getcwd(img_dir, sizeof(img_dir)) == NULL){
		printf("can't define current directory!\n");
		return EXIT_FAILURE;
	}
	strcat(img_dir,"/fc_imgs");
	mkdir(img_dir);
	memset(sys_com,0,strlen(sys_com));	
	while ( focused == 0 ){
		////////////////////////////////////////////////////
		sprintf(img_nme,"%s/FC%ld.fits",img_dir,time(NULL));	
		sprintf(sys_com, "camera takeimage %s %.2f 1 %d",img_nme,exposure,binning);
		mvsend_cmd2ccd(sys_com);
		sleep((int)exposure +2);
		memset(sys_com,0,strlen(sys_com));		
		sprintf(sys_com,"fitscopy '%s[%i:%i,%i:%i]' %s_subframe",img_nme,512-radius,512+radius,512-radius,512+radius,img_nme);
		system(sys_com);
		strcat(img_nme,"_subframe");
		/////////////////////////////////////////////////////
		fits_open_file(&f_in, img_nme, READONLY, &status_in);
		fits_get_img_param(f_in, 2, &bitpix, &naxis, naxes, &status_in);
		img_arr = (float *)malloc(naxes[0]*naxes[1]*sizeof(float)); 
		for ( fpx[1]=1; fpx[1] <= naxes[1]; fpx[1]++){
			curr_ptr = img_arr + (fpx[1]-1)*naxes[0];
			fits_read_pix(f_in, TFLOAT, fpx, naxes[0], NULL, curr_ptr, NULL, &status_in);
			if (arr_min == 0.0) arr_min = *(curr_ptr+0);
			while (fpx[0] < naxes[0]){
				if (*(curr_ptr+fpx[0]) >= arr_max) {arr_max = *(curr_ptr+fpx[0]);}
				if (*(curr_ptr+fpx[0]) > 0 && *(curr_ptr+fpx[0]) <= arr_min) arr_min = *(curr_ptr+fpx[0]);			
				fpx[0]+=1;
			}fpx[0] = 0; }
		printf("====================\nimage:%s\n",img_nme);
		printf("image_min:%.1f ; image_max:%.1f\n",arr_min,arr_max);
		///////////////////////////////////////////////////////////////////
		if (status_in != 0){ fits_report_error(stderr, status_in); break; } 
		else fits_close_file(f_in, &status_in);		
		if (arr_max >= ccd_limit_val*0.90){
			printf("values ​​exceeding the limit were found! Trying with exposure-=1");			
			exposure-=1; continue;		
		}
		for (fpx[1]=0;fpx[1]<naxes[1];fpx[1]++) for (fpx[0]=0;fpx[0]<naxes[0];fpx[0]++)
			*(img_arr+fpx[0]+naxes[0]*fpx[1]) = sqrt((*(img_arr+fpx[0]+naxes[0]*fpx[1])-arr_min)/(arr_max-arr_min));	
		gauss_filtr(img_arr,naxes,3);
		subtrct_MBG(img_arr,naxes[0]*naxes[1]);
		/////////////////////////////////////////////////////////	
		bn_mask = (float *)malloc(naxes[0]*naxes[1]*sizeof(float));
		printf("current threshold value: %.05f\n",(threshold = otsuThreshold(img_arr,bn_mask,naxes[0]*naxes[1]))); 		
		///////////////////////////////////////////		
		blob_el *blob_container = NULL, target_blob;
		blob_count = findallBlobs(bn_mask,naxes, &blob_container);
		if (blob_count == 0) { printf("no visible blobs on image"); break;}		
		else printf("found blobs count:%d\n",blob_count);		
		int targetb_num=0;
		float curr_r,min_r=naxes[0];
		if (bn_mask !=NULL) free(bn_mask);
		blob_el iter_blob;
		for (fpx[1]=0;fpx[1]<blob_count;fpx[1]++){
			iter_blob = blob_container[fpx[1]];
			curr_r = sqrt(pow((naxes[0]/2-iter_blob.x_cr),2)+pow((naxes[1]/2-iter_blob.y_cr),2));		
			if (curr_r <= min_r){
				min_r=curr_r; targetb_num=fpx[1];
		}	}
		target_blob = blob_container[targetb_num];		
		if (blob_container != NULL) free(blob_container);
		printf("target parameters:\n\tcenter:(%.02f;%.02f)\n\tsize(pxls):%i\n\tflux:%.02f\n",
			target_blob.x_cr, target_blob.y_cr, target_blob.countp, target_blob.flux); 
		/////////////////////////////
		x_cr = (int)target_blob.x_cr; 
		y_cr = (int)target_blob.y_cr;				
		ea_area = extract_EAA(img_arr,naxes,x_cr,y_cr,0.90);
		if (ea_area == -1){ printf("can't find effective ea_area!\n"); break; }
		target_row = ea_area*2+1;
		target_EAA = (float*)malloc(pow(target_row,2)*sizeof(float));
		x_cr = x_cr-ea_area; 
		y_cr = y_cr-ea_area;
		for (fpx[1]=0;fpx[1]<target_row;fpx[1]++) for(fpx[0]=0;fpx[0]<target_row;fpx[0]++){
			index = fpx[0]+fpx[1]*target_row;				
			target_EAA[index] = *(img_arr+(x_cr+fpx[0])+naxes[0]*(y_cr+fpx[1]));		
		}
		free(img_arr);
		//////////////////////
		naxes[0] = target_row;
		naxes[1] = naxes[0];
		strcat(img_nme,"_target");
		fits_create_file(&f_out, img_nme, &status_out);
		fits_create_img(f_out, FLOAT_IMG, 2, naxes, &status_out);
		fits_write_img(f_out, TFLOAT, 1, pow(target_row,2), target_EAA, &status_out);						
		if (status_out != 0){ fits_report_error(stderr, status_out); break; } 
		else fits_close_file(f_out, &status_out);
		///////////////////////////////////////
		HFD[0] = compute_HFD(target_EAA,naxes);
		printf("============\n HFD=%.04f\n============\n",HFD[0]);
		free(target_EAA);
		memset(img_nme,0,strlen(img_nme));
		focused = 1;				
		if (HFD[0] < HFD[1]){ 
			HFD[1] = HFD[0]; 
			printf("continue focusing...\n");			
		}else{	
			steps_per_request*=(-1);	
			focused = 1;
		}
		move_to = p_tat_info->fli_info.focuser_curr_position + steps_per_request;  
		sprintf(fli_cmd,"FOCUS %ld",move_to);
		mvsend_cmd2fli(fli_cmd);
		usleep((int)((float)steps_per_request*focuser_speed*1000000*2));
		printf("current_focuser_position:%ld\n",p_tat_info->fli_info.focuser_curr_position);
	}	     
	printf("=========\n FOCUSED \n=========\n");	
	remove_tat_info_shm( &shmid, p_tat_info);
	rmdir(img_dir);
	return EXIT_SUCCESS;
}
int equal_histo(float *image, long pxl_count){
	int x,y;
	float image_copy[pxl_count];
	for (x = 0; x < pxl_count; x++) image_copy[x] = image[x];
	histo *img_histo = histogram_maker(image_copy, pxl_count);
	for (x = 0; x < pxl_count; x++) for (y = 0; y < img_histo->hist_len; y++)
		if ( image[x] == img_histo->hist_val[y] ){ 
			image[x] += ((1.0-image[x])*img_histo->hist_cdf[y]);
			if (image[x] > 1.0) image[x] = 1.0;
			if (image[x] < 0.0) image[x] = 0.0;		
		} 
	free_histo(img_histo);
	return EXIT_SUCCESS;
}
float compute_HFD(float *target,long *naxes){
	float inRadius = naxes[0]/2, sum =0, sumDist =0, temp_val;
	int x,y, x_cr,y_cr;
	x_cr = naxes[0]/2; y_cr = x_cr;
	for (x=0;x<naxes[1];x++) for (y=0;y<naxes[0];y++)
		if (pow(x-x_cr,2) + pow(y-y_cr,2) <= pow(inRadius,2)){
			temp_val = *(target+x+naxes[0]*y); 
			sum+= temp_val;
			sumDist += temp_val * sqrt(pow((float)x-(float)x_cr, 2.0f) + pow((float)y-(float)y_cr, 2.0f));	
	}
	printf("sum(Vi*di)::%.02f sum(Vi)::%.02f\n",sumDist,sum);
	return (sum ? 2.0*sumDist/sum : sqrt(2.0)*inRadius);
}
float otsuThreshold(float *image, float *bn_mask, long pxl_count){
	int i,j, thresh, class1_pxlcount = 0;
	float best_sigma = 0.0, curr_intensicy_prob, second_class_prob, first_class_mean, second_class_mean, mean_delta, sigma;	
	float class1_Isum= 0.0, intensitySum=0.0, best_thresh = 0.0;
	float image_copy[pxl_count];
	for (i=0;i<pxl_count;i++) { image_copy[i] = image[i]; intensitySum+=image[i];}
	histo *hist = histogram_maker(image_copy, pxl_count);	
	for ( thresh = 0; thresh < hist->hist_len; ++thresh){
		class1_pxlcount += hist->hist_cnt[thresh];
		class1_Isum += hist->hist_val[thresh] * hist->hist_cnt[thresh];
		curr_intensicy_prob = class1_pxlcount / (float) pxl_count; 
		second_class_prob = 1.0 - curr_intensicy_prob;
		first_class_mean = class1_Isum / (float) class1_pxlcount;
		second_class_mean = (intensitySum - class1_Isum) / (float) (pxl_count - class1_pxlcount);
		mean_delta = first_class_mean - second_class_mean;
		sigma = curr_intensicy_prob * second_class_prob * mean_delta * mean_delta;
		if (sigma > best_sigma){
   			best_sigma = sigma;
    			best_thresh = hist->hist_val[thresh];
	}	}
	for (i = 0; i < pxl_count; i++) 
		bn_mask[i] = (image[i] >= best_thresh) ? (float)image[i] : 0;
	free_histo(hist);	
	return best_thresh;
}
histo * histogram_maker(float *image_copy, long pxl_count){
	int i,j,curr_prb;
	histo *new_histo = NULL;
	new_histo = (histo*)malloc(1*sizeof(histo));
	new_histo->hist_cnt = (int  *)malloc(1*sizeof(int));
	new_histo->hist_val = (float*)malloc(1*sizeof(float));
	new_histo->hist_cnt[0] = 1;
	new_histo->hist_val[0] = image_copy[0];
	qsort(image_copy, pxl_count, sizeof(float), comp);
	for (i=1,j=0;i<pxl_count;i++){
		if (image_copy[i] == image_copy[i-1]){ 
			new_histo->hist_cnt[j]+=1; continue; 
		} else {
			j+=1; 
			new_histo->hist_cnt = (int  *)realloc(new_histo->hist_cnt, (j+1)*sizeof(int));
			new_histo->hist_val = (float*)realloc(new_histo->hist_val, (j+1)*sizeof(float));
			new_histo->hist_cnt[j] = 1;
			new_histo->hist_val[j] = image_copy[i];
	}	}
	new_histo->hist_len = j+1;
	new_histo->hist_cdf = (float*)malloc(new_histo->hist_len*sizeof(float));
	for (i = 0; i < new_histo->hist_len; i++){
		new_histo->hist_cdf[i] = 0; 		
		for (j = 0; j <= i; j++) new_histo->hist_cdf[i] += (float)new_histo->hist_cnt[j]/pxl_count;	
	}
	return new_histo;
}
int free_histo(histo *curr_histo){
	free(curr_histo->hist_cnt);
	free(curr_histo->hist_val);
	free(curr_histo->hist_cdf);
	free(curr_histo);
	return EXIT_SUCCESS;
}
int gauss_filtr(float *img_arr, long *naxes, int kern_size){
	int x,y,i,j, kern_rad = kern_size/2.0;
	long new_naxes[2] = { naxes[0]+2*kern_rad, naxes[1]+2*kern_rad };
	float kernel[kern_size-1][kern_size-1], sigma=kern_rad/2.0, accumulation=0.0,conv_sum=0.0;
	float image_copy[new_naxes[0]*new_naxes[1]];
	memset( image_copy, 0, sizeof(float)*(new_naxes[0]*new_naxes[1]));
	for (y=0;y<kern_size;y++) for (x=0;x<kern_size;x++){
		kernel[y][x] = exp(-0.5*(pow((y-kern_rad)/sigma,2.0)+pow((x-kern_rad)/sigma,2.0)))/(2*M_PI*pow(sigma,2));
		accumulation+=kernel[y][x];
	}
	for (y=0;y<kern_size;y++) for (x=0;x<kern_size;x++) kernel[y][x]/=accumulation;
	for (y=0;y< naxes[0];y++) for (x=0;x< naxes[1];x++){
		if (x<kern_rad) 
			image_copy[((kern_rad-1)-x)+(y+kern_rad)*new_naxes[1]] = img_arr[x+y*naxes[1]];	
		if (naxes[0]-x<=kern_rad) 
			image_copy[(x+kern_rad-1)+2*(naxes[0]-x)+(y+kern_rad)*new_naxes[1]] = img_arr[x+y*naxes[1]];
		if (naxes[1]-y<=kern_rad) 
			image_copy[(x+kern_rad)+((y+kern_rad-1)+2*(naxes[1]-y))*new_naxes[1]] = img_arr[x+y*naxes[1]];
		if (y<kern_rad) 
			image_copy[(x+kern_rad)+((kern_rad-1)-y)*new_naxes[1]] = img_arr[x+y*naxes[1]];
		image_copy[(x+kern_rad)+(y+kern_rad)*new_naxes[1]] = img_arr[x+y*naxes[1]];
	}
	for ( y = kern_rad; y < naxes[0]; y++ ) for ( x = kern_rad; x < naxes[1]; x++ ) {
  		accumulation = 0.0;
  		conv_sum = 0.0; 
  		for (i=-kern_rad;i<=kern_rad;i++) for (j=-kern_rad;j<=kern_rad;j++) {
      			accumulation += image_copy[(x+j)+(y+i)*new_naxes[1]]*kernel[kern_rad+i][kern_rad+j];
      			conv_sum += kernel[kern_rad+i][kern_rad+j];
    		}
  	 	image_copy[x+y*new_naxes[1]] = accumulation/conv_sum;
	}
	for (y=0;y< naxes[0];y++) for (x=0;x< naxes[1];x++) 
		img_arr[x+y*naxes[1]] = image_copy[(x+kern_rad)+(y+kern_rad)*new_naxes[1]]; 
	return EXIT_SUCCESS; 		  
}
int extract_EAA(float *image,long *naxes, int x_cr, int y_cr, float accuracy){
	int edge_length=0,x,y;	
	double curr_I, prev_I;
	curr_I = *(image+x_cr+y_cr*naxes[0]);
	do{
		edge_length+=1;
		if (edge_length >= 15) { edge_length=-1; break; }
		prev_I = curr_I; curr_I=0;
		for (x=x_cr-edge_length;x<x_cr+edge_length;x++) for (y=y_cr-edge_length;y<y_cr+edge_length;y++) 
			curr_I+=*(image+x+naxes[0]*y);
		printf("prev_I/curr_I with R=%2d ==> %.02f%\n",edge_length, prev_I/curr_I);				
	} while (prev_I/curr_I < accuracy);
	return edge_length;
}
int findallBlobs(float *mask, long *naxes, blob_el **blob_array){
	int x,y,blob_count = 0, target_num=0;
	blob_el curr_blob={0,0,0};
	for (x=0;x<naxes[1];x++){		
		for (y=0;y<naxes[0];y++){
			curr_blob = (const blob_el){ 0 };
			if (blobDetector(mask,&curr_blob,naxes,x,y,1) == EXIT_SUCCESS) blob_count+=1;
			else continue;						
			if (blob_count == 1) 
				*blob_array = (blob_el*)malloc(sizeof(blob_el));
			else if (blob_count >= 2) 
				*blob_array = (blob_el*)realloc(*blob_array, blob_count * sizeof(blob_el));
			else continue;
			memcpy(*blob_array+blob_count-1,&curr_blob,sizeof(blob_el));
	}	}
	if (blob_count == 0) return EXIT_FAILURE;
	return blob_count;		
}
int blobDetector(float *mask, blob_el *blob, long *naxes, int x, int y, int flag){
	if (x < 0 || x >= naxes[0] || y < 0 || y >= naxes[1]) return EXIT_FAILURE;
	float curr_value = *(mask+x+naxes[0]*y);
	if (curr_value == 0) return EXIT_FAILURE;
	else{
		blob->x_cr+=x*curr_value;
		blob->y_cr+=y*curr_value;
		blob->flux+=curr_value;
		blob->countp+=1;
		*(mask+x+naxes[0]*y ) = 0;
		int i,conn[] = {1,0,-1,0,0,1,0,-1,1,1,-1,1,1,-1,-1,-1};
		for(i=0;i < (sizeof(conn)/2)/sizeof(int); i+=2){		
			blobDetector(mask,blob,naxes,x+conn[i],y+conn[i+1],0);
	}	}
	if (flag == 1){
		blob->x_cr=((*blob).x_cr/(*blob).flux);
		blob->y_cr=((*blob).y_cr/(*blob).flux);
	}
	return EXIT_SUCCESS;	
}
int subtrct_MBG(float *image, long pxl_count){ 
	int i;	
	float background, median, mean =0.0, image_copy[pxl_count];
	memcpy(image_copy, image, sizeof(float)*pxl_count);
	qsort(image_copy, pxl_count, sizeof(float),  comp);
	median = ( pxl_count%2 == 0 ) ? ((image_copy[pxl_count/2] + image_copy[pxl_count/2 - 1]) / 2.0) : image_copy[pxl_count/2];
	for (i =0; i< pxl_count; i++) mean += image_copy[i]; 
	mean/=pxl_count;
	background = 2.5*median - 1.5*mean;
	for (i =0; i< pxl_count; i++){
		image[i]-= background;
		if (image[i] < 0) image[i] = 0;
	}
	printf("mean:%f median:%f bg:%f\n",mean,median, background);
	return EXIT_SUCCESS;
}

