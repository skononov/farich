#ifndef RichSpectrum_h
# define RichSpectrum_h 1

# include <map>
# include <cmath> //for NAN value

class RichSpectrumIter;
class G4MaterialPropertyVector;

class RichSpectrum
{
	friend class RichSpectrumIter;
public:
	typedef std::map<double,double> data_type;
	typedef std::pair<double,double> point_type;
	typedef std::map<double,double>::iterator data_iterator;
	typedef std::map<double,double>::const_iterator data_const_iterator;

private:
	data_type data;
public:
	RichSpectrum() {}
	RichSpectrum(const char* datafile) { ReadFile(datafile); }
	RichSpectrum(double x1,double x2,double val) { data[x1]=val; data[x2]=val; }

	~RichSpectrum() { data.clear(); }

	int Size() const { return data.size(); }
	int IsEmpty() const { return data.empty(); }

	void Clear() { data.clear(); }

	bool AddEntry(double x,double val) { return data.insert(point_type(x,val)).second; }

	double& operator[](double x) { return (data.insert(point_type(x,0.0))).first->second; }

	void Set(double x1,double x2,double val) { Clear(); data[x1]=val; data[x2]=val; }

	int ReadFile(const char* datafile);

	void Scale(double c);

	double Evaluate(double x);

	void GetRange(double& x1,double& x2) {
		data_iterator first=data.begin(), last=data.end();
		x1=first->first;
		x2=(--last)->first;
	}

	void ReduceToRange(double x1,double x2,double defval=NAN);

	G4MaterialPropertyVector* GetMPV() const;

	void Print();
};

class RichSpectrumIter {
private:
	RichSpectrum::data_iterator iter;
	RichSpectrum* pspectrum;
public:
	RichSpectrumIter(RichSpectrum& spectrum) {
		pspectrum=&spectrum;
		iter=spectrum.data.begin();
	}
	void Init(RichSpectrum& spectrum) {
		pspectrum=&spectrum;
		iter=spectrum.data.begin();
	}
	void Reset() { iter=pspectrum->data.begin(); }
	RichSpectrum::point_type* operator()() {
		if(iter!=pspectrum->data.end())
			return (RichSpectrum::point_type*)&(*iter++);
		else
			return 0;
	}
};

#endif
