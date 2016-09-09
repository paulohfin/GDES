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

int contador, contador2;
uint32_t pto;
double cf[2];

typedef struct{
	double c[3];
}No;
std::vector<No> lista_centro;

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
double nvAngulo(double x, double y, double x1, double y1){
	double seno, coss;
	seno = asin((y1 - y) / sqrt(pow(x - x1, 2) + pow(y - y1, 2))) * 180 / 3.14159265;
	coss = acos((x1 - x) / sqrt(pow(x - x1, 2) + pow(y - y1, 2))) * 180 / 3.14159265;
	if(seno > 0) return coss;
	else return 360 - coss;
}
double ajuste(double beta2, double beta1, int var){
	if(((int)(beta2 - beta1) % 360 < var && (int)(beta2 - beta1) % 360 > -var) || ((int)(beta2 - beta1) % 360 > 360 - var && (int)(beta2 - beta1) % 360 < -360 + var)) return beta2;
	else return beta1;
}
double distancia2(double x1, double y1, double x2, double y2){
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}
double distancia(double x, double y, int T){
	double dist = sqrt(pow(x - cf[0], 2) + pow(y - cf[1], 2));
	if(dist < x || dist < y || dist < T - x || dist < T - y)
		return dist;
	else if(x < y && x < T - x && x < T - y){
		cf[0] = 0;
		cf[1] = y;
		return x;
	}
	else if(T - x < y && T - x < T - y){
		cf[0] = T;
		cf[1] = y;
		return T - x;
	}
	else if(y < T - y){
		cf[0] = x;
		cf[1] = 0;
		return y;
	}
	cf[0] = x;
	cf[1] = T;
	return T - y;
}
void nearest(double x, double y, ISpatialIndex* tree){
	double c[2];
	c[0] = x;
	c[1] = y;
	SpatialIndex::Point p = Point(&c[0], 2);
	MyVisitor vis;
	tree->nearestNeighborQuery(1, p, vis);
}
void Rua(ISpatialIndex** tree, FILE *fp, double x, double y, int r, int beta1, double x0, double y0, int T, int num){
	if(contador < -1 || contador2 < 0) return;

	int i, n;
	double beta2, c[2];
	No nv;

	std::vector<No> lista_aux;

	nearest(x, y, (*tree));
	nv.c[0] = x;
	nv.c[1] = y;
	nv.c[2] = beta1;
	beta2 = beta1;

	lista_aux.push_back(nv);
	while((nv.c[0] == x && nv.c[1] == y) || distancia(nv.c[0], nv.c[1], T) > r){		
		
		beta2 = ajuste(beta2 + (rand() % 20) - 10, beta1, 120);
		nv.c[0] += r * cos(3.1415926535 * beta2 / 180);
		nv.c[1] += r * sin(3.1415926535 * beta2 / 180);
		nv.c[2] = beta2;

		lista_aux.push_back(nv);
		nearest(nv.c[0], nv.c[1], (*tree));

		if(nv.c[0] < 0 || nv.c[0] > T || nv.c[1] < 0 || nv.c[1] > T){
			cf[0] = nv.c[0];
			cf[1] = nv.c[1];
			break;
		}	
	}
	lista_aux.pop_back();
	nv.c[0] = cf[0];
	nv.c[1] = cf[1];
	nv.c[2] = beta2;
	lista_aux.push_back(nv);

	if(lista_aux.size() < 3) return;

	contador--;
	contador2--;
	SpatialIndex::Point p;
	for(i = 0; i < lista_aux.size(); i++){
		c[0] = lista_aux[i].c[0];
		c[1] = lista_aux[i].c[1];

		p = Point(&c[0], 2);

		std::ostringstream os;
		os << p;
		std::string data = os.str();
		(*tree)->insertData(data.size() + 1, reinterpret_cast<const byte*>(data.c_str()), p, pto++);
	}
	
	if(contador > -1)		
		fprintf(fp,"\t\t{ \"type\": \"Feature\", \"properties\": { \"id\": %d }, \"geometry\": { \"type\": \"LineString\", \"coordinates\": [",contador);
	else
		fprintf(fp,"\t\t{ \"type\": \"Feature\", \"properties\": { \"id\": %d }, \"geometry\": { \"type\": \"LineString\", \"coordinates\": [",contador + 1);
		
	fprintf(fp," [ %lf , %lf ] ",lista_aux[0].c[0], lista_aux[0].c[1]);
	for(i = 1; i < lista_aux.size(); i++)
		fprintf(fp,", [ %lf , %lf ] ", lista_aux[i].c[0], lista_aux[i].c[1]);

	if(contador > -1)
		fprintf(fp," ] } }, \n");
	else
		fprintf(fp," ] } } \n");
	contador--;
	for(i = 0; i < lista_aux.size() && num > 1; i++){
		n = rand() % 6;
		if(n == 0 || n == 2)
			Rua(&(*tree), fp, lista_aux[i].c[0], lista_aux[i].c[1], r, lista_aux[i].c[2] - 90, x0, y0, T, num - 1);
		if(n == 1 || n == 2)
			Rua(&(*tree), fp, lista_aux[i].c[0], lista_aux[i].c[1], r, lista_aux[i].c[2] + 90, x0, y0, T, num - 1);
	}
}
int Avenida(ISpatialIndex** tree, FILE *fp, double x, double y, double r, int T, int R){
	double x1, y1, beta1, beta2, raio, c[3];
	int i, cont;
	No nv;
	if(x < 0) x1 = 0;
	else if(x > T) x1 = T;
	else x1 = x;
	if(y < 0) y1 = 0;
	else if(y > T) y1 = T;
	else y1 = y;
	std::vector<No> lista_aux;
	std::vector<No> lista_aux2;

	nearest(x1, y1, (*tree));

	cont = (int)distancia(x1, y1, T) / r;
	nv.c[0] = x1;
	nv.c[1] = y1;
	nv.c[2] = beta1;

	lista_aux.push_back(nv);
	beta1 = rand() % 360;
	beta2 = beta1;

	while((nv.c[0] == x1 && nv.c[1] == y1) || distancia(nv.c[0], nv.c[1], T) > r){			
		beta2 = ajuste(beta2 + (rand() % 60) - 30, beta1, 80);
		nv.c[0] += r * cos(3.1415926535 * beta2 / 180);
		nv.c[1] += r * sin(3.1415926535 * beta2 / 180);
		nv.c[2] = beta2;

		lista_aux.push_back(nv);
		cont--;
		if(cont <= 0){
			nearest(nv.c[0], nv.c[1], (*tree));
			cont = (int)distancia(nv.c[0], nv.c[1], T) / r;
		}
		if(distancia2(nv.c[0], nv.c[1], x, y) > R + rand() % (2 * R) && (nv.c[0] < 0 || nv.c[0] > T || nv.c[1] < 0 || nv.c[1] > T)){
			cf[0] = nv.c[0];
			cf[1] = nv.c[1];
			break;
		}
	}
	if(lista_aux.size() < 10) return 0;
	lista_aux.pop_back();
	nv.c[0] = cf[0];
	nv.c[1] = cf[1];
	nv.c[2] = beta2;
	lista_aux.push_back(nv);

	nv.c[0] = x1;
	nv.c[1] = y1;
	nv.c[2] = beta1 + 180;

	lista_aux2.push_back(nv);
	beta1 = rand() % 360;
	beta2 = beta1;

	while((nv.c[0] == x1 && nv.c[1] == y1) || distancia(nv.c[0], nv.c[1], T) > r){			
		beta2 = ajuste(beta2 + (rand() % 60) - 30, beta1, 80);
		nv.c[0] += r * cos(3.1415926535 * beta2 / 180);
		nv.c[1] += r * sin(3.1415926535 * beta2 / 180);
		nv.c[2] = beta2;

		lista_aux2.push_back(nv);
		cont--;
		if(cont <= 0){
			nearest(nv.c[0], nv.c[1], (*tree));
			cont = (int)distancia(nv.c[0], nv.c[1], T) / r;
		}
		if(distancia2(nv.c[0], nv.c[1], x, y) > R + rand() % (2 * R) && (nv.c[0] < 0 || nv.c[0] > T || nv.c[1] < 0 || nv.c[1] > T)){
			cf[0] = nv.c[0];
			cf[1] = nv.c[1];
			break;
		}
	}
	lista_aux2.pop_back();
	nv.c[0] = cf[0];
	nv.c[1] = cf[1];
	nv.c[2] = beta2;
	lista_aux2.push_back(nv);


	SpatialIndex::Point p;
	for(i = 0; i < lista_aux.size(); i++){
		c[0] = lista_aux[i].c[0];
		c[1] = lista_aux[i].c[1];

		p = Point(&c[0], 2);

		std::ostringstream os;
		os << p;
		std::string data = os.str();
		(*tree)->insertData(data.size() + 1, reinterpret_cast<const byte*>(data.c_str()), p, pto++);
	}
	for(i = 0; i < lista_aux2.size(); i++){
		c[0] = lista_aux2[i].c[0];
		c[1] = lista_aux2[i].c[1];

		p = Point(&c[0], 2);

		std::ostringstream os;
		os << p;
		std::string data = os.str();
		(*tree)->insertData(data.size() + 1, reinterpret_cast<const byte*>(data.c_str()), p, pto++);
	}
	
	if(contador > -1)		
		fprintf(fp,"\t\t{ \"type\": \"Feature\", \"properties\": { \"id\": %d }, \"geometry\": { \"type\": \"LineString\", \"coordinates\": [",contador);
	else
		fprintf(fp,"\t\t{ \"type\": \"Feature\", \"properties\": { \"id\": %d }, \"geometry\": { \"type\": \"LineString\", \"coordinates\": [",contador + 1);
	fprintf(fp," [ %lf , %lf ] ",lista_aux2[lista_aux2.size() - 2].c[0], lista_aux2[lista_aux2.size() - 2].c[1]);
	for(i = lista_aux2.size() - 2; i > 0; i--)
		fprintf(fp,", [ %lf , %lf ] ", lista_aux2[i].c[0], lista_aux2[i].c[1]);
	for(i = 0; i < lista_aux.size(); i++)
		fprintf(fp,", [ %lf , %lf ] ", lista_aux[i].c[0], lista_aux[i].c[1]);

	if(contador > -1)
		fprintf(fp," ] } }, \n");
	else
		fprintf(fp," ] } } \n");
	contador--;
	return 1;
}
void Rodovia(ISpatialIndex** tree, FILE *fp, double x, double y, double beta1, int T, int r){
	if(x < 0 || x > T || y < 0 || y > T) return;
	int i, raio, raio2, cont;
	double beta2, dist, c[2];

	std::vector<No> lista_aux;
	std::vector<No> lista_aux2;
	No nv;

	nearest(x, y, (*tree));
	dist = distancia(x, y, T);
	cont = (int)dist / r;
	nv.c[0] = x;
	nv.c[1] = y;
	beta2 = beta1;

	lista_aux.push_back(nv);

	while((nv.c[0] == x && nv.c[1] == y) || distancia(nv.c[0], nv.c[1], T) >= 1.5 * r){
		beta2 = ajuste(beta2 + (rand() % 30) - 15, beta1, 60);
		nv.c[0] += r * cos(3.1415926535 * beta2 / 180);
		nv.c[1] += r * sin(3.1415926535 * beta2 / 180);
		lista_aux.push_back(nv);
		cont--;
		if(cont <= 0){
			nearest(nv.c[0], nv.c[1], (*tree));
			dist = distancia(nv.c[0], nv.c[1], T);
			cont = (int)dist / r;
		}
		if(nv.c[0] < 0 || nv.c[0] > T || nv.c[1] < 0 || nv.c[1] > T){
			cf[0] = nv.c[0];
			cf[1] = nv.c[1];
			break;
		}
	}
	nv.c[0] = cf[0];
	nv.c[1] = cf[1];
	lista_aux.pop_back();
	lista_aux.push_back(nv);

	nv.c[0] = x;
	nv.c[1] = y;
	beta1 += 180;
	beta2 = beta1;
	cont = (int)dist / r;
	while((nv.c[0] == x && nv.c[1] == y) || distancia(nv.c[0], nv.c[1], T) >= 1.5 * r){
		beta2 = ajuste(beta2 + (rand() % 30) - 15, beta1, 45);
		nv.c[0] += r * cos(3.1415926535 * beta2 / 180);
		nv.c[1] += r * sin(3.1415926535 * beta2 / 180);
		lista_aux2.push_back(nv);
		cont--;
		if(cont <= 0){
			nearest(nv.c[0], nv.c[1], (*tree));
			dist = distancia(nv.c[0], nv.c[1], T);
			cont = (int)dist / r;
		}
		if(nv.c[0] < 0 || nv.c[0] > T || nv.c[1] < 0 || nv.c[1] > T){
			cf[0] = nv.c[0];
			cf[1] = nv.c[1];
			break;
		}
	}
	nv.c[0] = cf[0];
	nv.c[1] = cf[1];
	lista_aux2.pop_back();
	lista_aux2.push_back(nv);

	SpatialIndex::Point p;
	for(i = 0; i < lista_aux.size(); i++){
		c[0] = lista_aux[i].c[0];
		c[1] = lista_aux[i].c[1];

		p = Point(&c[0], 2);

		std::ostringstream os;
		os << p;
		std::string data = os.str();
		(*tree)->insertData(data.size() + 1, reinterpret_cast<const byte*>(data.c_str()), p, pto++);
	}
	for(i = 0; i < lista_aux2.size(); i++){
		c[0] = lista_aux2[i].c[0];
		c[1] = lista_aux2[i].c[1];

		p = Point(&c[0], 2);

		std::ostringstream os;
		os << p;
		std::string data = os.str();
		(*tree)->insertData(data.size() + 1, reinterpret_cast<const byte*>(data.c_str()), p, pto++);
	}
		
	if(contador > -1)		
		fprintf(fp,"\t\t{ \"type\": \"Feature\", \"properties\": { \"id\": %d }, \"geometry\": { \"type\": \"LineString\", \"coordinates\": [",contador);
	else
		fprintf(fp,"\t\t{ \"type\": \"Feature\", \"properties\": { \"id\": %d }, \"geometry\": { \"type\": \"LineString\", \"coordinates\": [",contador + 1);

	fprintf(fp," [ %lf , %lf ] ",lista_aux2[lista_aux2.size() - 1].c[0], lista_aux2[lista_aux2.size() - 1].c[1]);
	for(i = lista_aux2.size() - 2; i >= 0; i--)
		fprintf(fp,", [ %lf , %lf ] ", lista_aux2[i].c[0], lista_aux2[i].c[1]);
	for(i = 0; i < lista_aux.size(); i++)
		fprintf(fp,", [ %lf , %lf ] ", lista_aux[i].c[0], lista_aux[i].c[1]);	
	if(contador > -1)
		fprintf(fp," ] } }, \n");
	else
		fprintf(fp," ] } } \n");
	contador--;	
}
int main(int argc, char** argv){
	srand((unsigned)time(NULL));
	int T, r, r2, n, i, j, cont, sum, beta;
	double x, y, sigma;

	FILE *fp;
	char nome[50];
	strcpy(nome, argv[1]);
	strcat(nome, ".geojson");
	fp = fopen(nome,"w");
	T = atoi(argv[2]);
	contador = atoi(argv[3]);
	n = atoi(argv[4]);
	sum = 0;
	for(i = 0; i < n; i++){
		No nv;
		nv.c[0] = atoi(argv[5 + 3 * i]);
		nv.c[1] = atoi(argv[6 + 3 * i]);
		nv.c[2] = atoi(argv[7 + 3 * i]);
		sum += atoi(argv[7 + 3 * i]);
		lista_centro.push_back(nv);
	}
	sigma = atof(argv[5 + 3 * i]);
	pto = 1;
	cf[0] = rand() % T;
	cf[1] = rand() % T;

	r = 30;

	fprintf(fp,"{\n\t\"type\": \"FeatureCollection\",\n\t\"crs\": { \"type\": \"name\", \"properties\": { \"name\": \"urn:ogc:def:crs:OGC:1.3:CRS84\"} },\n\t\"features\": [\n");
	try{
		std::string baseName = "file";
		IStorageManager* diskfile = StorageManager::createNewDiskStorageManager(baseName, 4096);

		StorageManager::IBuffer* file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);

		id_type indexIdentifier;
		ISpatialIndex* tree = RTree::createNewRTree(*file, 0.7, 10, 10, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
		i = 0;
		cont = sqrt(sqrt(contador))/2;
		while(cont > 0){
			for(j = 0; j < lista_centro[i].c[2] && cont > 0; j++)
				Rodovia(&tree, fp, lista_centro[i].c[0] + rand() % (int)(sqrt(sqrt(T))), lista_centro[i].c[1] + rand() % (int)(sqrt(sqrt(T))), rand() % 360, T, r);
			i = (i + 1) % n;
			cont --;
		}
		cont = sqrt(sqrt(contador));
		for(i = 0; cont-- > 0; i++){
			i %= n;
			for(j = 0; j < lista_centro[i].c[2] && cont > 0; j++){
				beta = rand() % 360;
				r2 = sqrt(T) + rand() % (int)(T * lista_centro[i].c[2] / sum);
				Avenida(&tree, fp, lista_centro[i].c[0] + r2 * cos(3.1415926535 * beta / 180), lista_centro[i].c[1] + r2 * sin(3.1415926535 * beta / 180), r / 3, T, T * lista_centro[i].c[2] / sum);
			}
		}
		cont = contador;
		for(i = 0; contador > -1; i++){
			i %= n;
			contador2 = cont / sum;
			for(j = 0; j < lista_centro[i].c[2]; j++){
				j %= (int)lista_centro[i].c[2];
				r2 = (sqrt(T) +  rand() % T/n) * pow(1.0 * (rand() % T) / T, sigma * 10);;
				beta = rand() % 360;
				Rua(&tree, fp, lista_centro[i].c[0] + r2 * cos(beta * 3.1415926535 / 180), lista_centro[i].c[1] + r2 * sin(beta * 3.1415926535 / 180), r / 6, beta, lista_centro[i].c[0], lista_centro[i].c[1], T, lista_centro[i].c[2]);
			}
		}
		fprintf(fp,"\t]\n}");
		fclose(fp);
		delete tree;
		delete file;
		delete diskfile;
	}catch (Tools::Exception& e){
		std::cerr << "******ERROR******" << std::endl;
		std::string s = e.what();
		std::cerr << s << std::endl;
		return -1;
	}
	
	return 0;
}
