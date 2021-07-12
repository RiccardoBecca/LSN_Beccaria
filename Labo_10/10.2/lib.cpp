#include "lib.h"


int PBC(int a, int n_cities){
    if(a<n_cities) return a;
    else return a-n_cities;
}


Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

Cities::Cities(){
    vector<vector<double>> v(m_numero, vector<double>(2));
    
    for(int i = 0; i < m_numero; i++){
        
        //THIS IS FOR CITIES ON A CIRCLE
        //double theta= ((double) rand() / (RAND_MAX))*2*M_PI;
        //v[i][0] = cos(theta);
        //v[i][1] = sin(theta);
        
        //THIS IS FOR CITIES ON A SQUARE
        double x= ((double) rand() / (RAND_MAX));
        double y= ((double) rand() / (RAND_MAX));
        v[i][0] = x;
        v[i][1] = y;
    }
    m_city=v;
}

Cities::Cities(int square_circle, int numero, Random rand){
    
    ifstream input("Cities.txt");
    
    vector<vector<double>> v(numero, vector<double>(2));
    
    /*for(int i = 0; i < numero; i++){
        if(square_circle==0){
            //double theta= ((double) rand() / (RAND_MAX))*2*M_PI;
            double theta=rand.Rannyu(0,2*M_PI);
            v[i][0] = cos(theta);
            v[i][1] = sin(theta);
        }
        else{
            //double x= ((double) rand() / (RAND_MAX));
            //double y= ((double) rand() / (RAND_MAX));
            double x= rand.Rannyu();
            double y= rand.Rannyu();
            //double x= rand.Gauss(0,1);
            //double y= rand.Gauss(0,1);
            v[i][0] = x;
            v[i][1] = y;
        }
    }
     */
    
    for(int i=0; i< numero; i++){
        double x,y;
        input >> x >> y;
        v[i][0]=x;
        v[i][1]=y;
    }
    
    m_city=v; 
    m_numero=numero;
}


Cities::Cities(vector<vector<double>> v){
    m_city=v;
}

void Cities:: Print(void){
    for(int i=0; i< m_city.size(); i++){
        for(int j=0; j<m_city[i].size(); j++){
            cout << m_city[i][j] << " ";
        }
        cout << endl;
    }
}

Figlio::Figlio(int num){
    for(int i=0; i<num;i++){
        m_son.push_back(i+1);
    }
    random_shuffle ( m_son.begin()+1, m_son.end() );
    m_len=0;
}


 
void Figlio:: Print(void){
    for(int j=0; j<m_son.size();j++){
        cout << m_son[j] << " ";
    }
    cout << " Length = " << m_len << endl;
}

void Figlio::Copy(Figlio oldfiglio){
    for(int i=0; i<oldfiglio.Size(); i++){
        m_son[i]=oldfiglio.Get(i);
    }
    m_len=oldfiglio.Getlen();

}


Generazione::Generazione(int lunghezza, int num_figli, double cross, double m1, double m2, double m3, double m4, Random rand){
    for(int i=0; i<num_figli; i++){
        Figlio son(lunghezza);
        m_gen.push_back(son);
    }
    m_lunghezza=lunghezza;
    m_num_figlio= num_figli;
    m_prob_cross=cross;
    m_prob_m1=m1;
    m_prob_m2=m2;
    m_prob_m3=m3;
    m_prob_m4=m4;
    m_rand=rand;
}

void Generazione:: Print(void){
    for(int i=0; i<m_gen.size(); i++){
        m_gen[i].Print();
    }
}

void Figlio::Check(int numero_citta){
    int count=0;
    for(int i=0; i<numero_citta; i++){
        for(int j=0; j<numero_citta; j++){
            if(m_son[i]==m_son[j]){
                count+=1;
            }
        }
        count-=1;
    }
    if(count==0) cout << "This son is OK" << endl;
    else cout << "This son is NOT OK " << endl;
}

void Figlio::Measure_len(Cities city){
    double len=0;
    for(int i=0; i< m_son.size(); i++){
        int a=m_son[PBC(i,city.Getlen())];
        int b=m_son[PBC(i+1,city.Getlen())];
        len+=sqrt(pow(city.Get(a-1,0)-city.Get(b-1,0),2) + pow(city.Get(a-1,1)-city.Get(b-1,1),2)  );
    }
    m_len=len;
}

void Generazione::Check(void){
    for(int i=0; i<m_gen.size(); i++){
        m_gen[i].Check(m_lunghezza);
    }
}


void Generazione::Measure(Cities city){
    for(int i =0; i< m_gen.size(); i++){
        double len=0;
        for(int j=0; j<m_lunghezza; j++){
            //cout << "HET2 " << m_gen[i].Get(0) <<endl;
            //m_gen[i].Print();
            int a=m_gen[i].Get(PBC(j,m_lunghezza));
            //cout << "HET2" <<endl;
            int b=m_gen[i].Get(PBC(j+1,m_lunghezza));
            len+=sqrt(pow(city.Get(a-1,0)-city.Get(b-1,0),2) + pow(city.Get(a-1,1)-city.Get(b-1,1),2)  );
        }
        m_gen[i].Setlen(len);
    }
}

void Generazione::Sorting(void){
    sort(m_gen.begin(), m_gen.end(), sortcol);
}

bool sortcol(Figlio & v1, Figlio & v2 ) {
    return v1.Getlen() < v2.Getlen();
}



void Generazione::Crossing(Cities city, Generazione mygen){

    for(int i=0; i<m_num_figlio; i+=2){
        //int p= int(m_num_figlio*pow(((double) rand() / (RAND_MAX)), 2));
        //int q= int(m_num_figlio*pow(((double) rand() / (RAND_MAX)), 2));
        //double alpha=(double) rand() / (RAND_MAX);
        int p= int(m_num_figlio*pow((m_rand.Rannyu()), 2));
        int q= int(m_num_figlio*pow((m_rand.Rannyu()), 2));
        double alpha=m_rand.Rannyu();
        if(alpha<0.6){

            //int cut_cross= int(((double) rand() / (RAND_MAX))*m_lunghezza);
            int cut_cross= int((m_rand.Rannyu())*m_lunghezza);
            Figlio f1(m_lunghezza);
            Figlio f2(m_lunghezza);
    
            for(int t=0; t<m_lunghezza; t++){
                f1.Set(t, mygen.Get_son(p).Get(t));
            }
            for(int t=cut_cross; t<m_lunghezza; t++){
                bool written=false;
                for(int f=0; f<m_lunghezza; f++){
                    bool check=true;
                    for(int j=0; j<t; j++){
                        if(mygen.Get_son(q).Get(f)==f1.Get(j)){
                            check=false;
                        }
                    }
                    if(check==true & written==false){
                        f1.Set(t, mygen.Get_son(q).Get(f));
                        written=true;
                    }
                }
            }
            double len1=0;
            for(int j=0; j<m_lunghezza; j++){
                int a=f1.Get(PBC(j,m_lunghezza));
                int b=f1.Get(PBC(j+1,m_lunghezza));
                len1+=sqrt(pow(city.Get(a-1,0)-city.Get(b-1,0),2) + pow(city.Get(a-1,1)-city.Get(b-1,1),2)  );
            }
            f1.Setlen(len1);
    
            for(int t=0; t<m_lunghezza; t++){
                f2.Set(t, mygen.Get_son(q).Get(t));
            }
            for(int t=cut_cross; t<m_lunghezza; t++){
                bool written=false;
                for(int f=0; f<m_lunghezza; f++){
                    bool check=true;
                    for(int j=0; j<t; j++){
                        if(mygen.Get_son(p).Get(f)==f2.Get(j)){
                            check=false;
                        }
                    }
                    if(check==true & written==false){
                        f2.Set(t,mygen.Get_son(p).Get(f));
                        written=true;
                    }
                }
            }
            double len2=0;
            for(int j=0; j<m_lunghezza; j++){
                int a=f2.Get(PBC(j,m_lunghezza));
                int b=f2.Get(PBC(j+1,m_lunghezza));
                len2+=sqrt(pow(city.Get(a-1,0)-city.Get(b-1,0),2) + pow(city.Get(a-1,1)-city.Get(b-1,1),2)  );
            }
            f2.Setlen(len2);
    
            //if(f1.Getlen()+f2.Getlen() < m_gen[p].Getlen()+m_gen[q].Getlen()){
            if(0<1){
                for(int t=0; t<m_lunghezza; t++){
                    m_gen[i].Set(t, f1.Get(t));
                    m_gen[i+1].Set(t, f2.Get(t));
                }
                m_gen[i].Setlen(f1.Getlen());
                m_gen[i+1].Setlen(f2.Getlen());
            }
            else{
                for(int t=0; t<m_lunghezza;t++){
                    m_gen[i].Set(t, mygen.Get_son(p).Get(t));
                    m_gen[i+1].Set(t, mygen.Get_son(q).Get(t));
                }
                m_gen[i].Setlen(mygen.Get_son(p).Getlen());
                m_gen[i+1].Setlen(mygen.Get_son(q).Getlen());
            }
        }
    
    
        else{
            for(int t=0; t<m_lunghezza;t++){
                m_gen[i].Set(t, mygen.Get_son(p).Get(t));
                m_gen[i+1].Set(t, mygen.Get_son(q).Get(t));
            }
            m_gen[i].Setlen(mygen.Get_son(p).Getlen());
            m_gen[i+1].Setlen(mygen.Get_son(q).Getlen());
    
        }
    }
}
 
 
void Generazione::Mutation(Cities city){
    for(int i=0; i< m_gen.size(); i++){
        
        //double alpha=(double) rand() / (RAND_MAX);
        double alpha=m_rand.Rannyu();
        
        if(alpha<0.08){ //SWAP
            //int a = int(((double) rand() / (RAND_MAX))*(m_lunghezza-1) +1);
            //int b = int(((double) rand() / (RAND_MAX))*(m_lunghezza-1) +1);
            int a = int((m_rand.Rannyu())*(m_lunghezza-1) +1);
            int b = int((m_rand.Rannyu())*(m_lunghezza-1) +1);
            
            m_gen[i].Swap(a,b);
            m_gen[i].Measure_len(city);
        }
         
        //alpha=(double) rand() / (RAND_MAX);
        alpha=m_rand.Rannyu();
        if(alpha<0.08){ //PERMUTATION
            
            //int r = 1 + ( rand() % ( m_lunghezza-2 - 1 + 1 ) ); //numero compreso tra 1 e m_lunghezza -2
            //int n = 1 + ( rand() % ( int((m_lunghezza-r)*0.5) - 1 + 1 ) );
            //int m = 0 + ( rand() % ( m_lunghezza-2*n-r - 0 + 1 ) );
            int r=m_rand.Rannyu(1,m_lunghezza-2);
            int n=m_rand.Rannyu(1,(m_lunghezza-r)*0.5);
            int m=m_rand.Rannyu(0,m_lunghezza-2*n-r);
            Figlio copy(m_lunghezza);
            for(int j=0; j< m_gen[i].Size(); j++){
                copy.Set(j, m_gen[i].Get(j));
            }
            for(int j=0; j< r; j++){
                m_gen[i].Set(j, copy.Get(j));
            }
            for(int j =r; j<r+n;j++){
                m_gen[i].Set(j, copy.Get(j+n+m));
            }
            for(int j=r+n; j<r+n+m;j++){
                m_gen[i].Set(j, copy.Get(j));
            }
            for(int j=r+n+m;j<r+m+n+n;j++){
                m_gen[i].Set(j, copy.Get(j-n-m));
            }
            for(int j=r+m+n+n;j<m_gen[i].Size();j++){
                m_gen[i].Set(j, copy.Get(j));
            }
            m_gen[i].Measure_len(city);
        }
        
        //alpha=(double) rand() / (RAND_MAX);
        alpha=m_rand.Rannyu();
        if(alpha<0.08){ //SHIFT
            
            //int r = 1 + ( rand() % ( m_lunghezza-1 - 1 + 1 ) );
            //int n = 1 + ( rand() % ( m_lunghezza-r - 1 + 1 ) );
            //int m = 0 + ( rand() % ( m_lunghezza-n-r - 0 + 1 ) );
            int r=m_rand.Rannyu(1, m_lunghezza-1);
            int n=m_rand.Rannyu(1, m_lunghezza-r);
            int m=m_rand.Rannyu(0, m_lunghezza-n-r);
        
            Figlio copy(m_lunghezza);

            for(int j=0; j< m_gen[i].Size(); j++){
                copy.Set(j, m_gen[i].Get(j));
            }
            
            for(int j=0; j<r; j++){
                m_gen[i].Set(j, copy.Get(j));
            }
            for(int j=r; j<r+m; j++){
                m_gen[i].Set(j, copy.Get(j+n));
            }
            for(int j=r+m; j<r+n+m; j++){
                m_gen[i].Set(j, copy.Get(j-m));
            }
            for(int j=r+n+m; j<m_gen[i].Size();j++){
                m_gen[i].Set(j,copy.Get(j));
            }
            m_gen[i].Measure_len(city);
            
        }
        
        //alpha=(double) rand() / (RAND_MAX);
        alpha=m_rand.Rannyu();
        if(alpha<0.08){ //INVERSION
            Figlio copy(m_lunghezza);
            for(int j=0; j< m_gen[i].Size(); j++){
                copy.Set(j, m_gen[i].Get(j));
            }
            for(int j=1;j<m_lunghezza-1;j++){
                m_gen[i].Set(j, copy.Get(m_lunghezza-1-j));
            }
            m_gen[i].Measure_len(city);
        }
    }
}


void Generazione::Copy(Generazione oldgen){
    for(int i=0; i< m_gen.size(); i++){
        for(int j=0; j<oldgen.Get_son(i).Size(); j++){
            m_gen[i].Set(j, oldgen.Get_son(i).Get(j));
        }
        m_gen[i].Setlen(oldgen.Get_son(i).Getlen());
    }
}

double Generazione::Get_Best_Len(void){
    double min=m_gen[0].Getlen();
    for(int i=1;i<m_gen.size();i++){
        if(m_gen[i].Getlen()<min) min=m_gen[i].Getlen();
    }
    return min;
}
    
double Generazione::Get_Mean_Len(void){
    double mean=0;
    for(int i=0;i<int(0.5*m_gen.size());i++){
        mean+=m_gen[i].Getlen();
    }
    mean=mean/double(0.5*m_gen.size());
    return mean;
}

void Generazione::Change(Figlio f1, int pos){
    for(int i=0; i<f1.Size(); i++){
        m_gen[pos].Set(i, f1.Get(i));
    }
}
