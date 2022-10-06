class vars{ 
	 public:
		double G;
		double Q;
		std::complex<double> E;

		vars(){this->setTo(0.0);}

		void add(vars *v, vars *r){
			for(unsigned int i=0; i<sizeof(vars)/sizeof(double); i+=1){
				((double*)r)[i] = ((double*)this)[i] + ((double*)v)[i];
			}
		} 

		void mult(vars *r, double m){
			for(unsigned int i=0; i<sizeof(vars)/sizeof(double); i+=1){
				((double*)r)[i] = ((double*)this)[i]*m ;
			}
		}

		void setTo(double s){
			for(unsigned int i=0; i<sizeof(vars)/sizeof(double); i+=1){
				((double*)this)[i] = s;
			}
		}

};

enum vars_enum{
	G = 0,
	Q = 1,
	E = 2,
};
