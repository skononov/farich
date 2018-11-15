#ifndef RichAerogelMaterial_h
#define RichAerogelMaterial_h 1

# include <cmath>
# include <vector>
# include <unordered_map>
# include "globals.hh"
# include "G4MaterialPropertyVector.hh"
# include "RichParameters.hh"

class G4Material;

// Define injected hash object function for 
//  mapping refractive index to aerogel materials 
struct RichHashRI {
    size_t operator()(G4double x) const { return (size_t)rint(1e5*fabs(x)); }
};

struct RichEqualRI {
    constexpr bool operator()(const G4double &lhs, const G4double &rhs) const 
    {
        return rint(1e5*fabs(lhs)) == rint(1e5*fabs(rhs));
    }
};

class RichAerogelMaterial {
public:
	RichAerogelMaterial(G4double Lsc,const G4String& AbsorpDataFile,G4int chrom=kQuartzModel);

	~RichAerogelMaterial();

	void SetRefScatLength(G4double Lsc) { refScatLength=Lsc; initialized=false; }
	void SetAbsDataFile(G4String fn) { absorpDataFile=fn; initialized=false; }
	void SetChromaticity(G4int chrom) { chromaticity = chrom; }

	G4int Initialize(G4int verboseLevel=0);

	G4double GetScatLength(G4double PhotMom=kReferencePhotMom) const
	{ return rayleighData->GetProperty(PhotMom); }

	G4double GetAbsLength(G4double PhotMom=kReferencePhotMom) const
	{ return absorpData->GetProperty(PhotMom); }

	G4int GetChromaticity() const { return chromaticity; }

	G4Material* GetAerogelWithIndex(G4double ri,G4double pm=kReferencePhotMom);

private:
	G4double refScatLength;
	G4String absorpDataFile;
	G4int chromaticity;

	G4MaterialPropertyVector* absorpData;
	G4MaterialPropertyVector* rayleighData;
	G4bool initialized;

	typedef std::unordered_map< G4double, G4Material*, RichHashRI, RichEqualRI> mat_ri_map;
	mat_ri_map aerogels;
};

#endif

