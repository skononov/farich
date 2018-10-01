#ifndef Spectrum_h
# define Spectrum_h 1

# include <map>
# include <cmath> //for NAN value

class SpectrumIter;

class Spectrum
{
	friend class SpectrumIter;
public:
	typedef std::map<double,double> data_type;
	typedef std::pair<double,double> point_type;
	typedef std::map<double,double>::iterator data_iterator;
	typedef std::map<double,double>::const_iterator const_data_iterator;

private:
	data_type data;

public:
	Spectrum() {}
	Spectrum(const char* datafile) { ReadFile(datafile); }
	Spectrum(double x1,double x2,double val) { data[x1]=val; data[x2]=val; }
	Spectrum(Spectrum& sp) : data(sp.data) {}

	~Spectrum() { data.clear(); }

	int Size() const { return data.size(); }
	int IsEmpty() const { return data.empty(); }

	void Clear() { data.clear(); }

	bool AddEntry(double x,double val) { return data.insert(point_type(x,val)).second; }

	double& operator[](double x) { return (data.insert(point_type(x,0.0))).first->second; }

	void Set(double x1,double x2,double val) { Clear(); data[x1]=val; data[x2]=val; }

	int ReadFile(const char* datafile);

	void Scale(double c);

	double Evaluate(double x) const;

	void GetRange(double& x1,double& x2) const {
		const_data_iterator first=data.begin(), last=data.end();
		x1=first->first;
		x2=(--last)->first;
	}

	void ReduceToRange(double x1,double x2,double defval=NAN);

	void Print();
};

class SpectrumIter {
private:
	Spectrum::data_iterator iter;
	Spectrum* pspectrum;
public:
	SpectrumIter(Spectrum& spectrum) {
		pspectrum=&spectrum;
		iter=spectrum.data.begin();
	}
	void Init(Spectrum& spectrum) {
		pspectrum=&spectrum;
		iter=spectrum.data.begin();
	}
	void Reset() { iter=pspectrum->data.begin(); }
	Spectrum::point_type* operator()() {
		if(iter!=pspectrum->data.end())
			return (Spectrum::point_type*)&(*iter++);
		else
			return 0;
	}
};

#endif
