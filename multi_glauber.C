void multi_glauber()
{
	gSystem->Load("libMathMore");
	gSystem->CompileMacro("runglauber_v2.2_C", "k0");

	// const int NX = 7;
	// double x[NX] = {0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
	const int NX = 21;
	double x[NX] = {
		0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
		0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99
	};
	// const int NX = 1;
	// double x[NX] = {0.01};

	for (int ix = 0; ix < NX; ix++)
	{
		cout << endl;
		cout << "---------" << x[ix] << "---------" << endl;

		runAndSaveNtuple(100000, "p", "Au", 42, -1, 0.4, x[ix],
		                 Form("rootfiles/glauber_pau_snn42_beta141_x%03.0f_ntuple_100k.root", x[ix] * 100));

		runAndSaveNtuple(100000, "dh", "Au", 42, -1, 0.4, x[ix],
		                 Form("rootfiles/glauber_dau_snn42_beta141_x%03.0f_ntuple_100k.root", x[ix] * 100));

		runAndSaveNtuple(100000, "He3", "Au", 42, -1, 0.4, x[ix],
		                 Form("rootfiles/glauber_he3au_snn42_beta141_x%03.0f_ntuple_100k.root", x[ix] * 100));

		// runAndSaveNtuple(100000, "dh", "Au", 42, -1, 0.4, x[ix],
		//                  Form("rootfiles/glauber_dau_snn42_beta148_x%03.0f_ntuple_100k.root", x[ix]*100));
	}
}