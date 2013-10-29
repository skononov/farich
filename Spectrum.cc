#include "Spectrum.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using std::pair;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;

int Spectrum::ReadFile(const char* datafile)
{
	ifstream ifile(datafile);
	if (!ifile.is_open()) {
		cerr << datafile << ": No such file" << endl;
		return -1;
	}

	data.clear();

	double x, val;
	while (1) {
		ifile>>x>>val;
		if (ifile.fail() || ifile.eof()) break;
		pair<data_iterator,bool> result = data.insert(point_type(x,val));
		if (!result.second) {
			cerr<<datafile<<": entry ("<<x<<","<<val<<") duplicates point ("
				<<result.first->first<<","<<result.first->second<<")";
		}
	}
	ifile.close();

	return data.size();
}

void Spectrum::ReduceToRange(double x1,double x2,double defval)
{
	if (data.empty()) Set(x1,x2,defval);

	data_iterator first=data.begin(), last=data.end();
	last--;

	double xfirst=first->first, xlast=last->first;
	double yfirst=first->second, ylast=last->second;

	if (data.size()==1) //if only a single point defined set it for the given range
		Set(x1,x2,yfirst);

	if (x1>xlast || x2<xfirst) { //data are not defined for the given range
		Set(x1,x2,defval);
		return;
	}

	if (x1<xfirst) { //x1 below data range
		if( isnan(defval) )
			data[x1]=yfirst;
        else
			data[x1]=defval;
	} else { //loc1==0: x1 inside data range
		pair<data_iterator,bool> res=data.insert(point_type(x1,Evaluate(x1)));
		data.erase(data.begin(),res.first); //shrink data below x1
	}

	if (x2>xlast) { //x2 above data range
		if( isnan(defval) )
			data[x2]=ylast;
        else
			data[x2]=defval;
	} else { //loc2==0: x2 inside data range
		pair<data_iterator,bool> res=data.insert(point_type(x2,Evaluate(x2)));
		data.erase(++res.first,data.end()); //shrink data above x2
	}
}

void Spectrum::Scale(double c)
{
	if (data.empty()) return;

	data_iterator first=data.begin(), last=data.end();
	for(; first!=last; first++)
		first->second*=c;
}

double Spectrum::Evaluate(double x)
{
	if (data.empty()) return NAN;

	data_iterator first=data.begin(), last=data.end();
	last--;

	double val1=first->second, val2=last->second;

	data_iterator lower = data.lower_bound(x), upper=lower;

	if (lower==data.begin()) return val1; //falls below lower edge or on the first point

	if (lower==data.end()) return val2; //falls above upper edge

	if (lower->first==x) return lower->second; //exact match

	lower--; //take a point actually lower than x;

	double interval = upper->first-lower->first; //interval should be always greater than zero

	return (lower->second*(upper->first-x)+upper->second*(x-lower->first))/interval;
}

using std::setw;
using std::setprecision;
void Spectrum::Print()
{
	data_iterator first=data.begin(), last=data.end();
	for(; first!=last; first++)
		cout<<"  "<<setw(5)<<setprecision(4)<<first->first<<"  "<<first->second<<endl;
	cout.precision(0);
	cout.width(0);
}

