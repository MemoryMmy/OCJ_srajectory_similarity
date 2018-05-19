/* 
* date:   2017-07-29
* author: Mengyu Ma@National University of Defense Technology
* e-mail: mamengyu10@nudt.edu.cn
* description:
* 	This is a similar trajectories finding program we submitted to the ACM SIGSPATIAL Cup 2017. 
* In this competition, the goal is to find the similar trajectories of given trajectories by using 
* Fréchet distance as similarity measurement. In our method, we create spatial indexes for the first 
* and last points of the trajectories in the dataset.txt separately, which are used to filter out most 
* of the dissimilar trajectories and generate a much smaller candidate set. Then an Ordered-Coverage-Judge
* Fréchet distance algorithm is presented to search the similar trajectories form the candidate set.
*/
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <iostream> 
#include <stack>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/foreach.hpp>

#define MAX_PATH_SIZE 128
#define MIN_PATH_SIZE 8
#define MIN_ITEM_SIZE 8
#define MAX_NODE_SIZE 16



namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
namespace bgm = boost::geometry::model;


void Get_List(char* argv, char* result[], char* flag, int& count)
{
	char* string = strdup(argv);         		
	char* p;
	int i = 0;
	while((p = strsep(&string, flag)) != NULL)       
	{
		result[i] = p;
		i++;
	}
	result[i] = string;
	count = i;
}

bool Get_First2Item(char* argv, char* result[], char* flag)
{
	char* string = strdup(argv); 
	char* p;        		
	int i = 0;
	while((p = strsep(&string, flag))!= NULL&& i<2)       
	{
		result[i] = p;
		i++;
	}
	if (i<2) return false;
	return true;
}

void Get_Dir(char *pstring,char *pcrrentdir)
{
	int count=1;
    while((*(pcrrentdir++)=*(pstring++))) count++;
    while(*(--pcrrentdir)!='/' && count) count--;
    if(count)
		*pcrrentdir='\0';
	else
	{
		*(++pcrrentdir)='.';
		*(++pcrrentdir)='\0';
	}
}

void Get_Path(char *relatpath,char *dir,char *filepath)
{
	while((*(filepath++)=*(dir++)));
	*(filepath-1)='/';
	while((*(filepath++)=*(relatpath++)));
}

bool Make_Dir(const char *spathname)  
{  
    char dirname[256];  
    strcpy(dirname, spathname);  
    int i,len = strlen(dirname);  
    if(dirname[len-1]!='/')  
    strcat(dirname, "/");     
    len   =   strlen(dirname);    
    for(i=1;   i<len;   i++)  
    {  
        if(dirname[i]=='/')  
        {  
             dirname[i] = '\0';
             if(access(dirname,F_OK)!=0)    
                  if(mkdir(dirname,0755)<0)
					  return false;
             dirname[i] = '/';  
         }  
    }      
    return true;  
} 

int main(int argc, char* argv[])
{
	typedef bgm::d2::point_xy<double> point;
	typedef bgm::box<point> box;  
    typedef std::pair<point, unsigned> value;  

	char* datasetFile = argv[1];
    char* queryFile = argv[2];
	char* resultPath = argv[3];
	
		
	int myId, numProcs;
	FILE* file;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myId);
	MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
	
	//~ double t1,t2,t3,t4,t5,tx;
	//~ t1=MPI_Wtime();
	
	int datasetSize,querySize;
	char* datasetFileDir = new char[MAX_PATH_SIZE];
	char* queryFileDir = new char[MAX_PATH_SIZE];
	char* databuf;
	long len;
	int fd;

	fd=open(datasetFile,O_RDONLY);
	if(fd<0)
        printf("[ERROR] Can not open %s\n",datasetFile);  
	len = lseek(fd, 0, SEEK_END);  
	databuf = (char *) mmap( NULL,  len ,PROT_READ, MAP_PRIVATE, fd, 0 );
	char** datasetList= new char*[len/MIN_PATH_SIZE];
	Get_List(databuf,datasetList,(char*)"\n",datasetSize);
	munmap(databuf, len);
	close(fd);
	
	fd=open(queryFile,O_RDONLY);
	if(fd<0)
        printf("[ERROR] Can not open %s\n",queryFile);    
	len = lseek(fd, 0, SEEK_END);  
	databuf = (char *) mmap( NULL,  len ,PROT_READ, MAP_PRIVATE, fd, 0 );
	char** queryList= new char*[len/MIN_PATH_SIZE];
	Get_List(databuf,queryList,(char*)"\n",querySize);
	munmap(databuf, len);
	close(fd);

	Get_Dir(datasetFile,datasetFileDir);
	Get_Dir(queryFile,queryFileDir);
	
	if(myId==0)
		if(!Make_Dir(resultPath))  
			printf("[ERROR] Can not create dir %s\n",resultPath);

	//~ tx=MPI_Wtime();
	//~ printf("[TEST] Successfully. cores:%d coreid:%d cost:<total>%f</total> s\n",numProcs,myId, tx-t1);
	
	double** datasetPolyXList=new double*[datasetSize];
	double** datasetPolyYList=new double*[datasetSize];
	double* datasetPolyMaxxList=(double*)malloc(sizeof(double)*datasetSize);
	double* datasetPolyMinxList=(double*)malloc(sizeof(double)*datasetSize);
	double* datasetPolyMaxyList=(double*)malloc(sizeof(double)*datasetSize);
	double* datasetPolyMinyList=(double*)malloc(sizeof(double)*datasetSize);
	int* datasetPolyLList=(int*)malloc(sizeof(int)*datasetSize);
	
	char* pathbuf = new char[MAX_PATH_SIZE];
	
	//create rtree index
	bgi::rtree< value, bgi::quadratic<MAX_NODE_SIZE> > rtreeHead;
	bgi::rtree< value, bgi::quadratic<MAX_NODE_SIZE> > rtreeEnd;
	std::vector<value> headList;
	std::vector<value> endList;	

	int interval =datasetSize/numProcs;
    for(int itmp=interval*myId;itmp<datasetSize+interval*myId;itmp++)
    {
		int i=itmp%datasetSize;
		int trajSize;
		Get_Path(datasetList[i],datasetFileDir,pathbuf);
		fd=open(pathbuf,O_RDONLY);   
		len = lseek(fd, 0, SEEK_END);  
		if(fd<0)
		    printf("[WARNING] Can not open %s\n",pathbuf);
		if((int)len>0)
		{
			char** trajList= new char*[len/MIN_ITEM_SIZE];
			databuf = (char *) mmap( NULL,  len ,PROT_READ, MAP_PRIVATE, fd, 0 );	
			Get_List(databuf,trajList,(char*)"\n",trajSize);
			munmap(databuf, len);
			close(fd);
			
			datasetPolyXList[i]=(double*)malloc(sizeof(double)*(trajSize-1));
			datasetPolyYList[i]=(double*)malloc(sizeof(double)*(trajSize-1));
			datasetPolyLList[i]=0;
			bool isFirst=true;
			for(int j=1;j<trajSize;j++)
			{
				char* pxy[2];
				
				if(Get_First2Item(trajList[j],pxy,(char*)" "))
				{
					double px = atof(pxy[0]);
					double py = atof(pxy[1]);
					if(isFirst)
					{
						isFirst=false;
						datasetPolyMaxxList[i] = px;
						datasetPolyMinxList[i] = px;
						datasetPolyMaxyList[i] = py;
						datasetPolyMinyList[i] = py;
					}
					else
					{
						if(datasetPolyMaxxList[i] < px) datasetPolyMaxxList[i] = px;
						if(datasetPolyMinxList[i] > px) datasetPolyMinxList[i] = px;
						if(datasetPolyMaxyList[i] < py) datasetPolyMaxyList[i] = py;
						if(datasetPolyMinyList[i] > py) datasetPolyMinyList[i] = py;
						
					}
					datasetPolyXList[i][j-1] = px;
					datasetPolyYList[i][j-1] = py;
					datasetPolyLList[i]++;
				}
			
			}
			headList.push_back(std::make_pair(point(datasetPolyXList[i][0],datasetPolyYList[i][0]), i));
			endList.push_back(std::make_pair(point(datasetPolyXList[i][datasetPolyLList[i]-1],
			datasetPolyYList[i][datasetPolyLList[i]-1]), i));
		}   
	}
	//~ t2 = MPI_Wtime();
	//~ printf("[DATA LOAD] Successfully. cores:%d coreid:%d cost:<total>%f</total> s\n",numProcs,myId, t2-t1);
	rtreeHead.insert(headList.begin(),headList.end());
	rtreeEnd.insert(endList.begin(),endList.end());
	
	//~ t3 = MPI_Wtime();
	//~ printf("[INDEXES BUILD] Successfully. cores:%d coreid:%d cost:<total>%f</total> s\n",numProcs,myId, t3-t2);
    
    for(int i=myId;i<querySize;i+=numProcs)
    {
		int trajSize;
		char* querybuf[2];
		double FDBound;
		if(Get_First2Item(queryList[i],querybuf,(char*)" "))
		{
			Get_Path(querybuf[0],queryFileDir,pathbuf);
			FDBound=strtof(querybuf[1],NULL);
			fd=open(pathbuf,O_RDONLY);   
			len = lseek(fd, 0, SEEK_END);
			if(fd<0)
				printf("[WARNING] Can not open %s\n",pathbuf);  
			databuf = (char *) mmap( NULL,  len ,PROT_READ, MAP_PRIVATE, fd, 0 );
			char** trajList= new char*[len/MIN_ITEM_SIZE];
			Get_List(databuf,trajList,(char*)"\n",trajSize);
			munmap(databuf, len);
			close(fd);

			double* queryPolyX=(double*)malloc(sizeof(double)*(trajSize-1));
			double* queryPolyY=(double*)malloc(sizeof(double)*(trajSize-1));
			double queryPolyMaxx,queryPolyMinx,queryPolyMaxy,queryPolyMiny;
			int queryPolyL=0;
			bool isFirst=true;
			for(int j=1;j<trajSize;j++)
			{	
				char* pxy[2];
				
				if(Get_First2Item(trajList[j],pxy,(char*)" "))
				{
					double px = atof(pxy[0]);
					double py = atof(pxy[1]);
					if(isFirst)
					{
						isFirst=false;
						queryPolyMaxx = px;
						queryPolyMinx = px;
						queryPolyMaxy = py;
						queryPolyMiny = py;
					}
					else
					{
						if(queryPolyMaxx < px) queryPolyMaxx = px;
						if(queryPolyMinx > px) queryPolyMinx = px;
						if(queryPolyMaxy < py) queryPolyMaxy = py;
						if(queryPolyMiny > py) queryPolyMiny = py;
					}
					queryPolyX[j-1] = px;
					queryPolyY[j-1] = py;
					queryPolyL++;
				}
			}
			//query using head and end
			box headBuffer;
			box endBuffer;
			bg::buffer(bg::return_envelope<box>(point(queryPolyX[0],queryPolyY[0])),headBuffer,FDBound);
			bg::buffer(bg::return_envelope<box>(point(queryPolyX[queryPolyL-1],queryPolyY[queryPolyL-1])),endBuffer,FDBound);

			std::vector<value> headResult;
			std::vector<value> endResult;
			
			rtreeHead.query(bgi::intersects(headBuffer), std::back_inserter(headResult));
			rtreeEnd.query(bgi::intersects(endBuffer), std::back_inserter(endResult));
			int* indexList = new int[headResult.size()];
			int indexLen=0;
			BOOST_FOREACH(value const& vh, headResult)
				BOOST_FOREACH(value const& ve, endResult)
					if(vh.second==ve.second)
					{
						indexList[indexLen]=vh.second;
						indexLen++;
					}
			

			char *writebuf= new char[indexLen*40];
			writebuf[0]='\0';

			for(int j=0;j<indexLen;j++)
			{	

				int v=indexList[j];
				if((datasetPolyMaxxList[v]<=queryPolyMaxx+FDBound)
				&&(datasetPolyMinxList[v]>=queryPolyMinx-FDBound)
				&&(datasetPolyMaxyList[v]<=queryPolyMaxy+FDBound)
				&&(datasetPolyMinyList[v]>=queryPolyMiny-FDBound))
				{
					// Frechet Filter
					std::stack<int>  cIndexStack;
					std::stack<int>  pHeadStack;
					std::stack<int>  pEndStack;
					bool isSimilar=false;
					double squareBound = FDBound * FDBound;
					int pHead = 0;
					int pEnd = datasetPolyLList[v]-1;
					double dX=queryPolyX[0]-datasetPolyXList[v][0];
					double dY=queryPolyY[0]-datasetPolyYList[v][0];
					double squareDistance =dX*dX+dY*dY;
					if (squareDistance<=squareBound)//
					{
						dX=queryPolyX[queryPolyL-1]-datasetPolyXList[v][datasetPolyLList[v]-1];
						dY=queryPolyY[queryPolyL-1]-datasetPolyYList[v][datasetPolyLList[v]-1];
						squareDistance =dX*dX+dY*dY;
						if(squareDistance<=squareBound)
						{
							isSimilar=true;
							int matN = datasetPolyLList[v] * queryPolyL;
							bool (*flagMat)[queryPolyL] = (bool(*)[queryPolyL])malloc(matN);
							memset(flagMat,true,matN);
							
							for(int di=1;di<datasetPolyLList[v];di++)
							{
								dX=queryPolyX[0]-datasetPolyXList[v][di];
								dY=queryPolyY[0]-datasetPolyYList[v][di];
								squareDistance =dX*dX+dY*dY;
								if(squareDistance > squareBound)
								{
									pEnd=di;
									break;							
								}
							}
							
							
							for(int qi=1;qi<queryPolyL-1;qi++)
							{
								
								isSimilar = false;
								bool hasStrategy = false;
								if(flagMat[pHead][qi])
								{
									flagMat[pHead][qi]=false;
									int pEndCopy = pEnd;

									for(int di=pHead;di<=pEndCopy;di++)
									{
										dX=queryPolyX[qi]-datasetPolyXList[v][di];
										dY=queryPolyY[qi]-datasetPolyYList[v][di];
										squareDistance = dX*dX+dY*dY;
										if(squareDistance <= squareBound)
										{
											if(!isSimilar)
											{
												if(hasStrategy)
												{
													cIndexStack.push(qi);
													pHeadStack.push(pHead);								
													pEndStack.push(pEnd);
												}
												pHead = di;
												isSimilar = true;
											}
											if(di==pEndCopy)
											{
												int dj;
												for(dj=di+1;dj<datasetPolyLList[v];dj++)
												{
													dX=queryPolyX[qi]-datasetPolyXList[v][dj];
													dY=queryPolyY[qi]-datasetPolyYList[v][dj];
													squareDistance = dX*dX+dY*dY;
													if(squareDistance > squareBound)
													{
														
														pEnd=dj;
														break;							
													}
													
												}
												if(dj==datasetPolyLList[v])
												{
													pEnd=dj-1;
												}
												hasStrategy = true;
											}
										}
										else
										{
											if (isSimilar)
											{
												pEnd = di;
												hasStrategy = true;
												isSimilar = false;
											}
										}
									}
								}
								if(!hasStrategy)
								{
									if(cIndexStack.empty())
									{
										isSimilar=false;
										break;
									}
									else
									{
										qi = cIndexStack.top();
										pHead = pHeadStack.top();								
										pEnd = pEndStack.top();
										cIndexStack.pop();
										pHeadStack.pop();								
										pEndStack.pop();
									}
								}
								else
									isSimilar=true;
							}
							if(pEnd<datasetPolyLList[v]-1)
							{
								for(int di =pEnd;di<datasetPolyLList[v]-1;di++)
								{
									dX=queryPolyX[queryPolyL-1]-datasetPolyXList[v][di];
									dY=queryPolyY[queryPolyL-1]-datasetPolyYList[v][di];
									squareDistance = dX*dX+dY*dY;
									if(squareDistance > squareBound)
									{
										isSimilar=false;
										break;							
									}
								}
							}
							free(flagMat);
						}
							
					}
					if(isSimilar)
					{
						strcat(writebuf,datasetList[v]);
						strcat(writebuf,(char*)"\n");
					}
					
					
				}
			}
			
			char* resultFile = new char[MAX_PATH_SIZE];
			sprintf(resultFile,"%s/result-%05d.txt",resultPath,i);
			file=fopen(resultFile, "wt"); 
			fwrite(writebuf, 1, strlen(writebuf), file);
			fclose(file);
			
			free(queryPolyX);
			free(queryPolyY);
			delete[] writebuf;
			delete[] resultFile;
			delete[] indexList;
			                 
		}	
    }
	
	//~ MPI_Barrier(MPI_COMM_WORLD);
	//~ t4 = MPI_Wtime();
	//~ printf("[QUERY] Successfully. cores:%d coreid:%d cost:<total>%f</total> s\n",numProcs,myId, t4-t3);
	
	//~ t5 = MPI_Wtime();
	//~ printf("[TOTAL] Successfully. cores:%d coreid:%d cost:<total>%f</total> s\n",numProcs,myId, t5-t1);
	
	//~ if(myId==0)
	//~ {
		//~ printf("[RESULT] Successfully. cores:%d cost:<total>%f</total> s\n",numProcs, t5-t1);
	//~ }
		
	MPI_Finalize();
  
}
