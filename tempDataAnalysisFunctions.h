//This function gets the vector data[] and return the autocorrelation in function of
//separation


void extract_mean_sigma(int displ,double *data, int max, int size, double *meandeltav, double *sigma){
	int i,j;
	double temp;
	double *input = NULL;
	input = malloc( size*sizeof(double));

	for(j=0;j<max;j++){
		for(i=0;i<size;i++){
			input[i]=data[j*size+i];
			}
		
		for(i=0;i<(size-displ);i++){
			temp=input[i+displ]-input[i];
			*meandeltav+= temp;
			}

		*meandeltav = *meandeltav/((double)(size-displ)*(double)max);

		for(i=0;i<(size-displ);i++){
			temp=input[i+displ]-input[i];
			temp=temp-(*meandeltav);
			*sigma += pow(temp,2);
			}
	}

	*sigma=*sigma/((double)(size-displ)*(double)max);
	*sigma=sqrt(*sigma);

	free(input);
}

void histogram(FILE *fout,int displ, double *data, double *N, int size, int Nbins, double width, double *hflat){
	int max, aux, qflat = 4;
	int i, j;
	double temp, meandeltav,sum,sigma_aux,sigma=0.;
	double lastbin, firstbin, lbin;
	double auxflat1 = 0., auxflat2 = 0.;

	double *input = NULL;
	input = malloc( size*sizeof(double));
	float *count = NULL;
	count = malloc(Nbins*sizeof(float));

	for(i=0;i<Nbins;i++) {count[i]=0.;}
	*hflat = 0.;

	max = (int)( (*N)/size);

	extract_mean_sigma(displ, data, max, size, &meandeltav, &sigma);
	printf("meandeltav=%lf\nsigma=%lf\n",meandeltav,sigma);

	firstbin=-1.*width;
	lastbin=width;
	lbin=(lastbin-firstbin)/Nbins;


	for(j=0;j<max;j++){
		for(i=0;i<size;i++){
			input[i]=data[j*size+i];
			}
		
		for(i=0;i<(size-displ);i++){
			temp=(input[i+displ]-input[i]);
			//printf("temp = %lf\n",temp);
			auxflat1 += pow(fabs(temp-meandeltav),qflat);
			auxflat2+= pow(fabs(temp-meandeltav),3);
			temp=(temp-meandeltav)/sigma;
			aux=(int)((temp-firstbin)/lbin);
			//printf("aux = %d\n",aux);
			if (aux<0) {printf("aux negativo %d\n",aux);}
			if(aux>Nbins) {printf("Alerta de aux > Nbins!!\n");}
			if(aux==Nbins) {count[(int)(Nbins-1)]++;}
			else {count[aux]++;}
		}			
	}

	for(j=0;j<Nbins;j++){
		count[j]=count[j]/((double)(size-displ)*(double)max*lbin);
//		//sum+=count[j]*lbin;
		}		

	for(i=0;i<Nbins;i++) {fprintf(fout,"%lf\t%2.9lf\n",firstbin+i*lbin,count[i]);}
//	printf("sigma=%lf\n",(sigma));

	auxflat1 = auxflat1/((double)(size-displ)*(double)max);
	auxflat2 = auxflat2/((double)(size-displ)*(double)max);
	auxflat2 = pow(auxflat2,(double)qflat/3.);
	*hflat = auxflat1/auxflat2;

	free(input);
	free(count);
}
//--------------------------------------------------------------------
void extract_mean_sigma_from_file(int displ, FILE *fin, int max, int size, double *meandeltav, double *sigma, double *kurt){
	int i,j;
	double temp;
	float tempf;
	double *input = NULL;
	input = malloc( size*sizeof(double));

	fscanf(fin,"%f\n",&tempf);
	printf("temp = %f\n",tempf);

	for(j=0;j<max;j++){
		i=0;
		while(fscanf(fin,"%f\n",&tempf)!=EOF && i < size){
			input[i]=(double)tempf;
			i++;
			}
		
		for(i=0;i<(size-displ);i++){
			temp=input[i+displ]-input[i];
			*meandeltav+= temp;
			}

		*meandeltav = *meandeltav/((double)(size-displ)*(double)max);

		for(i=0;i<(size-displ);i++){
			temp=input[i+displ]-input[i];
			temp=temp-(*meandeltav);
			*sigma += pow(temp,2);
			*kurt += pow(temp,4);
			}
	}

	*sigma=*sigma/((double)(size-displ)*(double)max);
	*sigma=sqrt(*sigma);
	*kurt=*kurt/((double)(size-displ)*(double)max);
	*kurt = *kurt / pow(*sigma,4);

	free(input);
	rewind(fin);
}

//----------------------------------------------------------------------
void histogram_from_file(FILE *fout,int displ, FILE *fin, double *N, int size, int Nbins, double width, double *hflat, FILE *foutkurt){
	int max, aux, qflat = 4;
	int i, j;
	float tempf;
	double temp, meandeltav,sum,sigma_aux,sigma=0., kurt = 0.;
	double lastbin, firstbin, lbin;
	double auxflat1 = 0., auxflat2 = 0.;

	double *input = NULL;
	input = malloc( size*sizeof(double));
	double *count = NULL;
	count = malloc(Nbins*sizeof(double));

	for(i=0;i<Nbins;i++) {count[i]=0.;}
	*hflat = 0.;

	max = (int)( (*N)/size);
	printf("Ate aqui ok 1!\n");
	extract_mean_sigma_from_file(displ, fin, max, size, &meandeltav, &sigma, &kurt);
	printf("meandeltav=%lf\nsigma=%lf\nkurt=%lf\n",meandeltav,sigma,kurt);
	fprintf(foutkurt,"%d\t%lf\n",displ,kurt);

	firstbin=-1.*width;
	lastbin=width;
	lbin=(lastbin-firstbin)/Nbins;

	printf("Ate aqui ok 2!\n");
	for(j=0;j<max;j++){
		i=0;
		while(fscanf(fin,"%f\n",&tempf)!=EOF && i < size){
			input[i]=(double)tempf;			
			i++;
			//if(i==10) {printf("input[10] = %f\n",input[i]);}				
			}

		for(i=0;i<(size-displ);i++){
			temp=(input[i+displ]-input[i]);
			//printf("temp = %lf\n",temp);
			auxflat1 += pow(fabs(temp-meandeltav),qflat);
			auxflat2+= pow(fabs(temp-meandeltav),3);
			temp=(temp-meandeltav)/sigma;
			aux=(int)((temp-firstbin)/lbin);
			//printf("aux = %d\n",aux);
			if (aux<0) {printf("aux negativo %d\n",aux);}
			if(aux>Nbins) {printf("Alerta de aux > Nbins!!\n");}
			if(aux==Nbins) {count[(int)(Nbins-1)]++;}
			else {count[aux]++;}
		}			
	}
	for(j=0;j<Nbins;j++){
		count[j]=count[j]/((double)(size-displ)*(double)max*lbin);
//		//sum+=count[j]*lbin;
		}		
  
	for(i=0;i<Nbins;i++) {fprintf(fout,"%lf\t%2.9lf\n",firstbin+i*lbin,count[i]);}
//	printf("sigma=%lf\n",(sigma));

	auxflat1 = auxflat1/((double)(size-displ)*(double)max);
	auxflat2 = auxflat2/((double)(size-displ)*(double)max);
	auxflat2 = pow(auxflat2,(double)qflat/3.);
	*hflat = auxflat1/auxflat2;

	free(input);
	free(count);
	rewind(fin);
}

//------------------------------------------------------------------

void rearrange_y(double *data, double *data_rearr, double N){
	int i,j,k,l=0;
	int linesize = 1024;

	for(k=0;k<N/(linesize*linesize);k++){
		for(j=0; j< linesize; j++){
			for(i=0;i<linesize;i++){	
				data_rearr[l] = data [i*linesize + j + k*linesize*linesize];
				l++;
				}
			}
		}
	}


void structure_function(FILE *out, double *data, int displ, int size, double N, int order_min, int order_max){
	
	int i,j,k,max;
	double sum[order_max - order_min];
//	double *sum = NULL;
//	sum = malloc( (order_max - order_min) * sizeof(double));
	double temp;
	double *input = NULL;
	input = malloc(size*sizeof(double));

	for(k=0;k<(order_max - order_min);k++) {sum[k] = 0.;}

	printf("Hello! (1)\n");

	max = (int) ((N)/size);

	for(j=0;j<max;j++){
		for(i=0;i<size;i++){
			input[i]=data[j*size+i];
			}
		
		for(i=0;i<(size-displ);i++){
			temp=(input[i+displ]-input[i]);
			for(k=0;k<=(order_max - order_min);k++) {sum[k] += pow(temp,k+order_min);}
			}			
	}

	printf("Hello! (2)\n");

	for(k=0;k<=(order_max - order_min);k++) {sum[k] = sum[k] / ((double)(max*(size-displ)));}

	printf("Hello! (3)\n");

	fprintf(out,"%d",displ);
//	for(k=0;k<=(order_max - order_min);k++) {fprintf(out,"\t%lf",fabs(sum[k]));}
	for(k=0;k<=(order_max - order_min);k++) {fprintf(out,"\t%lf",sum[k]);}
	fprintf(out,"\n");

	printf("Hello! (4)\n");

	free(input);
//	free(sum);
	printf("Hello! (5)\n");

}

void structure_function_file(FILE *out, FILE *fdata, int displ, int size, double N, int order_min, int order_max){
	
	int i,j,k,max;
	double sum[order_max - order_min + 1];
//	double *sum = NULL;
//	sum = malloc( (order_max - order_min) * sizeof(double));
	double temp;
	double *input = NULL;
	input = malloc(size*sizeof(double));

	for(k=0;k<=(order_max - order_min);k++) {sum[k] =0.;}

	printf("Hello! (1)\n");

	max = (int) ((N)/size);
	rewind(fdata);

	for(j=0;j<max;j++){
		for(i=0;i<size;i++){
			fscanf(fdata,"%lf\n",&input[i]);
//			input[i]=data[j*size+i];
			}
		
		for(i=0;i<(size-displ);i++){
			temp=(input[i+displ]-input[i]);
			for(k=0;k<=(order_max - order_min);k++) {sum[k] += pow(temp,k+order_min);}
			}			
	}

	rewind(fdata);

	printf("Hello! (2)\n");

	for(k=0;k<=(order_max - order_min);k++) {sum[k] = sum[k] / ((double)(max*(size-displ)));}

	printf("Hello! (3)\n");

	fprintf(out,"%d",displ);
//	for(k=0;k<=(order_max - order_min);k++) {fprintf(out,"\t%lf",fabs(sum[k]));}
	for(k=0;k<=(order_max - order_min);k++) {fprintf(out,"\t%lf",sum[k]);}
	fprintf(out,"\n");

	printf("Hello! (4)\n");

	free(input);
//	free(sum);
	printf("Hello! (5)\n");

}

