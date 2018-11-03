void getqe(int x,int y)
{
	double wl, pde;

	wl=(x-88.)*500./785.+200.;
	pde=(y-478.)/(20.-478.)*40.;

	cout<<wl<<"    "<<pde<<endl;
}
