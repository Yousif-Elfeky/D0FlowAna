#ifndef calmean_hh
#define calmean_hh
// #include <math>
/***This class is used for calculate the mean value at each fill***/

class calmean {
  public:
    calmean();
    calmean(const calmean &cm);
    virtual ~calmean(){}
    void fill(float num);
    float mean() ;
    int count() {
      return mCount;
    }
    float sum() const {
      return mSum;
    }
    void setsum(float a);
    void setcount(int a);
    void set(float s, int i); //set sum and counts
    void print();
  private:
    float mSum;
    float mAvg;
    int mCount;
};

inline calmean::calmean() :mSum(0),mAvg(0),mCount(0) {
}
inline calmean::calmean(const calmean &cm){
  mSum = cm.mSum;
  mCount = cm.mCount;
  if (mCount!=0) mAvg=1.0*mSum/mCount;
  else mAvg = 0;
}
inline void calmean::fill(float num){
  mCount++;
  mSum+=num;
}
inline float calmean::mean() {
  if (mCount!=0) mAvg = mSum*1.0/mCount;
  return mAvg;
}
inline void calmean::print(){
  if (mCount!=0) mAvg=1.0*mSum/mCount;
  cout<<"sum: "<<mSum<<", mean: "<<mAvg<<", total counts: "<<mCount<<endl;
}
inline void calmean::set(float s, int i){
  mCount = i;
  mSum = s;
  if (mCount==0) {
    cout<<"Warn: Counts is "<<0<<endl;
    cout<<"Warn: Reset the sum to  "<<0<<endl;
    mSum = 0;
  }
}
inline void calmean::setsum(float a){
  if (mCount==0) {
    cout<<"Warn: Counts is "<<0<<endl;
    cout<<"Reset the sum to  "<<0<<endl;
    mSum = 0;
  }
  else mSum = a;
}
inline void calmean::setcount(int a){
  if (a<0) cout<<"Error! The Counts is negative!"<<endl;
  else if (a==0 && mCount!=0) cout<<"Warn: Counts=0 while the previous sum is not zero. Setting counts fail."<<endl;
  else mCount = a;
}

#endif
