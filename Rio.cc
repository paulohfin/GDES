#include <iostream>
#include <spatialindex/capi/sidx_api.h>
#include <spatialindex/capi/sidx_impl.h>
#include <spatialindex/capi/sidx_config.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

using namespace std;
using namespace SpatialIndex;

int contador, ver;
uint32_t pto;
double cf[2];

typedef struct{
	double c[2];
}No;

typedef struct{
	std::vector<No> vertices;
}coord;
std::vector<coord> lista_centro;

class MyVisitor : public IVisitor{

	public:
		MyVisitor(){ } 

		void visitNode(const INode& n) {}

		void visitData(const IData& d){
			int id = d.getIdentifier();
			SpatialIndex::IShape* shape;
			d.getShape(&shape);

			SpatialIndex::Point center;
			shape->getCenter(center);
			cf[0] = center.m_pCoords[0];
			cf[1] = center.m_pCoords[1];
		}

		void visitData(std::vector<const IData*>& v) {}
};
void verificar(double x, double y, int T){
	double dist = sqrt(pow(x - cf[0], 2) + pow(y - cf[1], 2));
	if(dist < x || dist < y || dist < T - x || dist < T - y)
		return;
	else if(x < y && x < T - x && x < T - y){
		cf[0] = 0;
		cf[1] = y;
		return;
	}
	else if(T - x < y && T - x < T - y){
		cf[0] = T;
		cf[1] = y;
		return;
	}
	else if(y < T - y){
		cf[0] = x;
		cf[1] = 0;
		return;
	}
	cf[0] = x;
	cf[1] = T;
}
double nvAngulo(double x, double y, double x1, double y1){
	double seno, coss;
	seno = asin((y1 - y) / sqrt(pow(x - x1, 2) + pow(y - y1, 2))) * 180 / 3.14159265;
	coss = acos((x1 - x) / sqrt(pow(x - x1, 2) + pow(y - y1, 2))) * 180 / 3.14159265;
	if(seno > 0) return coss;
	else return 360 - coss;
}
void nearest(double x, double y, ISpatialIndex* tree){
	double c[2];
	c[0] = x;
	c[1] = y;
	SpatialIndex::Point p = Point(&c[0], 2);
	MyVisitor vis;
	tree->nearestNeighborQuery(1, p, vis);
}
double ajuste(double x0, double y0, double x1, double y1, double beta){
	double beta1, beta2;
	beta1 = nvAngulo(x0, y0, x1, y1);
	beta2 = nvAngulo(x1, y1, cf[0], cf[1]);
	if(abs(beta - beta1) > 85 && abs(beta - beta1) < 285) return beta2;
	if(abs(beta - beta2) > 85 && abs(beta - beta2) < 285) return beta2;
	return beta;
}
void Agregar(ISpatialIndex** tree, FILE *fp, double x, double y, int r, int T){
	if(x < 0 || x >= T || y < 0 || y >= T) return; 
	int i;
	std::vector<No> lista_aux;
	double beta;
	No nv;

	nv.c[0] = x;
	nv.c[1] = y;

	lista_aux.push_back(nv);

	nearest(x, y, (*tree));
	beta = nvAngulo(x, y, cf[0], cf[1]);

	nv.c[0] += r * cos(3.1415926535 * beta / 180);
	nv.c[1] += r * sin(3.1415926535 * beta / 180);

	verificar(nv.c[0], nv.c[1], T);
	
	while(sqrt(pow(nv.c[0] - cf[0], 2) + pow(nv.c[1] - cf[1], 2)) > r){
		beta = ajuste(x, y, nv.c[0], nv.c[1], beta + (rand() % 60) - 30);
		nv.c[0] += r * cos(3.1415926535 * beta / 180);
		nv.c[1] += r * sin(3.1415926535 * beta / 180);
		lista_aux.push_back(nv);
	}
	
	if(lista_aux.size() > 10){
		double c[2];
		for(i = 0; i < lista_aux.size(); i++){
			c[0] = lista_aux[i].c[0];
			c[1] = lista_aux[i].c[1];

			SpatialIndex::Point p = Point(&c[0], 2);

			std::ostringstream os;
			os << p;
			std::string data = os.str();
			(*tree)->insertData(data.size() + 1, reinterpret_cast<const byte*>(data.c_str()), p, pto++);
		}
		
		if(contador > -1)
			fprintf(fp,"\t\t{ \"type\": \"Feature\", \"properties\": { \"id\": %d }, \"geometry\": { \"type\": \"LineString\", \"coordinates\": [",contador);
		else fprintf(fp,"\t\t{ \"type\": \"Feature\", \"properties\": { \"id\": %d }, \"geometry\": { \"type\": \"LineString\", \"coordinates\": [",contador + 1);
		
		fprintf(fp," [ %lf , %lf ] ",lista_aux[0].c[0], lista_aux[0].c[1]);
		for(i = 1; i < lista_aux.size(); i++)
			fprintf(fp,", [ %lf , %lf ] ", lista_aux[i].c[0], lista_aux[i].c[1]);
		fprintf(fp,", [ %lf , %lf ] ", cf[0], cf[1]);
		
		if(contador > -1)
			fprintf(fp," ] } }, \n");
		else
			fprintf(fp," ] } } \n");
		contador--;
	}
}
void Distribuir(ISpatialIndex** tree, FILE *fp, int n, int T, double sigma, int r){
	int i, r2;
	double x, y, beta;
	i = 0;
	while(i < lista_centro[n].vertices.size() && contador > -2){
		r2 = (rand() % ((int)sqrt(2) * T)) * pow(1.0 * (rand() % T) / T, sigma * 10);
		beta = 1.0 * (rand() % (360 * T)) / T;
		x = lista_centro[n].vertices[i].c[0] + r2 * cos(3.14159265 * beta / 180);
		y = lista_centro[n].vertices[i].c[1] + r2 * sin(3.14159265 * beta / 180);
		Agregar(&(*tree), fp, x, y, r, T);
		i++;
	}
}
int main(int argc, char** argv){
	srand((unsigned)time(NULL));
	int T, r, beta, n, i, j;
	double r2, sigma;
	FILE *fp;
	char nome[50];
	strcpy(nome, argv[1]);
	strcat(nome,".geojson");
	fp = fopen(nome,"w");
	T = atoi(argv[2]);
	contador = atoi(argv[3]);
	n = atoi(argv[4]);
	for(i = 0; i < n; i++){
		coord nv;
		No nv2;
		nv2.c[0] = atoi(argv[5 + 3 * i]);
		nv2.c[1] = atoi(argv[6 + 3 * i]);
		nv.vertices.push_back(nv2);
		lista_centro.push_back(nv);
	}
	sigma = atof(argv[5 + 3 * i]);
	r = 5;
	pto = 1;
	cf[0] = 0;
	cf[1] = 0;

	fprintf(fp,"{\n\t\"type\": \"FeatureCollection\",\n\t\"crs\": { \"type\": \"name\", \"properties\": { \"name\": \"urn:ogc:def:crs:OGC:1.3:CRS84\"} },\n\t\"features\": [\n");
	try{
		std::string baseName = "file";
		IStorageManager* diskfile = StorageManager::createNewDiskStorageManager(baseName, 4096);

		StorageManager::IBuffer* file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);

		id_type indexIdentifier;
		ISpatialIndex* tree = RTree::createNewRTree(*file, 0.7, 10, 10, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);

		j = 0;
		while(contador > -2){
			Distribuir(&tree, fp, j, T, sigma, r);
			j = (j + 1) % n;
		}

		delete tree;
		delete file;
		delete diskfile;
	}catch (Tools::Exception& e){
		std::cerr << "******ERROR******" << std::endl;
		std::string s = e.what();
		std::cerr << s << std::endl;
		return -1;
	}
	fprintf(fp,"\t]\n}");
	fclose(fp);
	return 0;
}
