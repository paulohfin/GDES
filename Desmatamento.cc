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

int contador;
uint32_t line;

typedef struct{
	double c[2];
}coord;

typedef struct{
	std::vector<coord> vertices;
}No;

std::vector<No> lista;
std::vector<No> lista_centro;
std::vector<int> lid;

class MyVisitor : public IVisitor{
	public:
		MyVisitor(){ } 

		void visitNode(const INode& n) {}

		void visitData(const IData& d){
			lid.push_back(d.getIdentifier() - 1);
		}

		void visitData(std::vector<const IData*>& v) {}
};
void nearest(double x, double y, ISpatialIndex* tree){
	double c[2];
	c[0] = x;
	c[1] = y;
	SpatialIndex::Point p = Point(&c[0], 2);
	MyVisitor vis;
	tree->nearestNeighborQuery(100, p, vis);
}
int verificar(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
	double m = ((y4-y2) * (x3-x4) * (x2-x1) + x2 * (y2 - y1) * (x3-x4) - x4*(y3-y4) * (x2-x1)) / ((x3-x4) * (y2-y1) - (y3-y4) * (x2-x1));
    	
	if(((m >= x1 && m <= x2) || (m <= x1 && m >= x2)) && ((m >= x3 && m <= x4) || (m <= x3 && m >= x4)))
		return 1;
	return 0;
}
int verificar2(double x, double y, int id, int raio){
	int sum = 0, i;
	for(i = 0; i < lista[id].vertices.size(); i++)
		if(sqrt((x - lista[id].vertices[i].c[0]) * (x - lista[id].vertices[i].c[0]) + (y - lista[id].vertices[i].c[i]) * (y - lista[id].vertices[i].c[i]))< raio)
			sum++;
	if(sum > 0) return 1;
	else return 0;
}
void Agregar(ISpatialIndex** tree, FILE *fp, double x, double y, int T, int r){
	if(x < 0 || x > T || y < 0 || y > T) return;

	int i, j, k, raio, cont, raio2, n;

	double beta, c[2];

	lid.clear();

	nearest(x, y, (*tree));
	No nv;
	coord nv2;
	raio = 2 + rand() % r;
	for(j = 0; j < lid.size(); j++){
		cont = 0;
		if(verificar2(x, y, lid[j], raio) == 1){
			for(i = 0; i < lista[lid[j]].vertices.size(); i++)
				cont += verificar(x, y, 0, 0, lista[lid[j]].vertices[i].c[0], lista[lid[j]].vertices[i].c[1], lista[lid[j]].vertices[(i + 1) % lista[lid[j]].vertices.size()].c[0], lista[lid[j]].vertices[(i + 1) % lista[lid[j]].vertices.size()].c[1]);
			if(cont % 2 == 1)
				return;
		}
	}

	for(beta = 0;beta < 360;){
		raio2 = 5 + rand() % raio;
		nv2.c[0] = x + raio2 * cos(3.1415926535 * beta / 180);
		nv2.c[1] = y + raio2 * sin(3.1415926535 * beta / 180);
		if(nv2.c[0] < 0) nv2.c[0] = 0;
		if(nv2.c[0] > T - 1) nv2.c[0] = T - 1;
		if(nv2.c[1] < 0) nv2.c[1] = 0;
		if(nv2.c[1] > T - 1) nv2.c[1] = T - 1;

		nv.vertices.push_back(nv2);
		beta += rand() % 120;
	}
	for(k = 0; k < lid.size(); k++)
		for(j = 0;j < lista[lid[k]].vertices.size(); j++)
			for(i = 0; i < nv.vertices.size(); i++)
				if(verificar(nv.vertices[(i + 1) % nv.vertices.size()].c[0], nv.vertices[(i + 1) % nv.vertices.size()].c[1], nv.vertices[i].c[0], nv.vertices[i].c[1], lista[lid[k]].vertices[(j + 1) % lista[lid[k]].vertices.size()].c[0], lista[lid[k]].vertices[(j + 1) % lista[lid[k]].vertices.size()].c[1], lista[lid[k]].vertices[j].c[0], lista[lid[k]].vertices[j].c[1]) == 1) return;

	lista.push_back(nv);
	SpatialIndex::Point p;
	c[0] = x;
	c[1] = y;

	p = Point(&c[0], 2);

	std::ostringstream os;
	os << p;
	std::string data = os.str();
	(*tree)->insertData(data.size() + 1, reinterpret_cast<const byte*>(data.c_str()), p, line++);
	
	if(contador > -1)
		fprintf(fp,"\t\t{ \"type\": \"Feature\", \"properties\": { \"id\": %d }, \"geometry\": { \"type\": \"LineString\", \"coordinates\": [",contador);
	else fprintf(fp,"\t\t{ \"type\": \"Feature\", \"properties\": { \"id\": %d }, \"geometry\": { \"type\": \"LineString\", \"coordinates\": [",contador + 1);
	
	fprintf(fp," [ %lf , %lf ] ", nv.vertices[0].c[0], nv.vertices[0].c[1]);
	for(i = 1; i < nv.vertices.size(); i++)
		fprintf(fp,", [ %lf , %lf ] ", nv.vertices[i].c[0], nv.vertices[i].c[1]);
	fprintf(fp," , [ %lf , %lf ] ", nv.vertices[0].c[0], nv.vertices[0].c[1]);
	
	if(contador > -1)
		fprintf(fp," ] } }, \n");
	else
		fprintf(fp," ] } } \n");
	contador--;
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
		Agregar(&(*tree), fp, x, y, T, r);
		i++;
	}
}
int main(int argc, char** argv){
	srand((unsigned)time(NULL));
	int T, r, n, i, j;
	double r2, sigma, beta;
	FILE *fp;

	char nome[50];
	strcpy(nome, argv[1]);
	strcat(nome,".geojson");
	fp = fopen(nome,"w");
	T = atoi(argv[2]);
	contador = atoi(argv[3]);
	r = atoi(argv[4]);
	n = atoi(argv[5]);
	for(i = 0; i < n; i++){
		No nv;
		coord nv2;
		j = 0;
		while(j < atoi(argv[8 + 3 * i])){
			r2 = (rand() % ((int)sqrt(2) * T));
			beta = 1.0 * (rand() % (360 * T)) / T;
			nv2.c[0] = atoi(argv[6 + 3 * i]) + r2 * cos(3.14159265 * beta / 180);
			nv2.c[1] = atoi(argv[7 + 3 * i]) + r2 * sin(3.14159265 * beta / 180);
			nv.vertices.push_back(nv2);
			j++;
		}
		lista_centro.push_back(nv);
	}
	sigma = atof(argv[6 + 3 * i]);
	line = 1;

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

