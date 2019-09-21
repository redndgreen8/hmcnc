

#include <iostream>
#include <time.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

/* Generates a sequence of states and emissions from a Markov model */
static void generate(const int &_N,const cv::Mat &_TRANS,const cv::Mat &_EMIS, const cv::Mat &_INIT, cv::Mat &seq, cv::Mat &states)
{
    seq = cv::Mat(1,_N,CV_32S);
    states = cv::Mat(1,_N,CV_32S);
    int n_states = _TRANS.rows;
    cv::Mat cumulative_emis(_EMIS.size(),CV_64F);
    for (int r=0;r<cumulative_emis.rows;r++)
        cumulative_emis.at<double>(r,0) = _EMIS.at<double>(r,0);
    for (int r=0;r<cumulative_emis.rows;r++)
        for (int c=1;c<cumulative_emis.cols;c++)
            cumulative_emis.at<double>(r,c) = cumulative_emis.at<double>(r,c-1) + _EMIS.at<double>(r,c);
    cv::Mat cumulative_trans(_TRANS.size(),CV_64F);
    for (int r=0;r<cumulative_trans.rows;r++)
        cumulative_trans.at<double>(r,0) = _TRANS.at<double>(r,0);
    for (int r=0;r<cumulative_trans.rows;r++)
        for (int c=1;c<cumulative_trans.cols;c++)
            cumulative_trans.at<double>(r,c) = cumulative_trans.at<double>(r,c-1) + _TRANS.at<double>(r,c);
    cv::Mat cumulative_init(_INIT.size(),CV_64F);
    cumulative_init.at<double>(0,0) = _INIT.at<double>(0,0);
    for (int c=1;c<cumulative_init.cols;c++)
        cumulative_init.at<double>(0,c) = cumulative_init.at<double>(0,c-1) + _INIT.at<double>(0,c);
    double r_init,r_trans,r_emis;
    r_init = (double) rand()/RAND_MAX;
    int last_state;
    for (int c=0;c<cumulative_init.cols;c++)
        if (r_init <= cumulative_init.at<double>(0,c))
        {
            last_state = c;
            break;
        }
    for (int t=0;t<_N;t++)
    {
        r_trans = (double)rand()/RAND_MAX;
        for (int i=0;i<cumulative_trans.cols;i++)
            if (r_trans <= cumulative_trans.at<double>(last_state,i))
            {
                states.at<int>(0,t) = i;
                break;
            }
        r_emis = (double)rand()/RAND_MAX;
        for (int i=0;i<cumulative_emis.cols;i++)
        {
            if (r_emis <= cumulative_emis.at<double>(states.at<int>(0,t),i))
            {
                seq.at<int>(0,t) = i;
                break;
            }
        }
        last_state = states.at<int>(0,t);
    }
}
static void generate(const int &_N,const int &_M, const cv::Mat &_TRANS,const cv::Mat &_EMIS, const cv::Mat &_INIT, cv::Mat &seq, cv::Mat &states)
	{
		seq = cv::Mat(_M,_N,CV_32S);
		states = cv::Mat(_M,_N,CV_32S);
		for (int i=0;i<_M;i++)
		{
			cv::Mat seq_,states_;
			generate(_N,_TRANS,_EMIS,_INIT,seq_,states_);
			for (int t=0;t<_N;t++)
			{
				seq.at<int>(i,t) = seq_.at<int>(0,t);
				states.at<int>(i,t) = states_.at<int>(0,t);
			}
		}

	}


void correctModel(cv::Mat &TRANS, cv::Mat &EMIS, cv::Mat &INIT)
{
    double eps = 1e-30;
    for (int i=0;i<EMIS.rows;i++)
        for (int j=0;j<EMIS.cols;j++)
            if (EMIS.at<double>(i,j)==0)
                EMIS.at<double>(i,j)=eps;
    for (int i=0;i<TRANS.rows;i++)
        for (int j=0;j<TRANS.cols;j++)
            if (TRANS.at<double>(i,j)==0)
                TRANS.at<double>(i,j)=eps;
    for (int i=0;i<INIT.cols;i++)
        if (INIT.at<double>(0,i)==0)
            INIT.at<double>(0,i)=eps;
    double sum;
    for (int i=0;i<TRANS.rows;i++)
    {
        sum = 0;
        for (int j=0;j<TRANS.cols;j++)
            sum+=TRANS.at<double>(i,j);
        for (int j=0;j<TRANS.cols;j++)
            TRANS.at<double>(i,j)/=sum;
    }
    for (int i=0;i<EMIS.rows;i++)
    {
        sum = 0;
        for (int j=0;j<EMIS.cols;j++)
            sum+=EMIS.at<double>(i,j);
        for (int j=0;j<EMIS.cols;j++)
            EMIS.at<double>(i,j)/=sum;
    }
    sum = 0;
    for (int j=0;j<INIT.cols;j++)
        sum+=INIT.at<double>(0,j);
    for (int j=0;j<INIT.cols;j++)
        INIT.at<double>(0,j)/=sum;
}
static void printModel(const cv::Mat &TRANS,const cv::Mat &EMIS,const cv::Mat &INIT)
{
    std::cout << "\nTRANS: \n";
    for (int r=0;r<TRANS.rows;r++)
    {
        for (int c=0;c<TRANS.cols;c++)
            std::cout << TRANS.at<double>(r,c) << " ";
        std::cout << "\n";
    }
    std::cout << "\nEMIS: \n";
    for (int r=0;r<EMIS.rows;r++)
    {
        for (int c=0;c<EMIS.cols;c++)
            std::cout << EMIS.at<double>(r,c) << " ";
        std::cout << "\n";
    }
    std::cout << "\nINIT: \n";
    for (int r=0;r<INIT.rows;r++)
    {
        for (int c=0;c<INIT.cols;c++)
            std::cout << INIT.at<double>(r,c) << " ";
        std::cout << "\n";
    }
    std::cout << "\n";
}


static void getUniformModel(const int &n_states,const int &n_observations, cv::Mat &TRANS,cv::Mat &EMIS,cv::Mat &INIT)
{
    TRANS = cv::Mat(n_states,n_states,CV_64F);
    TRANS = 1.0/n_states;


    INIT = cv::Mat(1,n_states,CV_64F);
    INIT = 1.0/n_states;

    EMIS = cv::Mat(n_states,n_observations,CV_64F);
    EMIS = 1.0/n_observations;
}

/* Calculates maximum likelihood estimates of transition and emission probabilities from a sequence of emissions */
void train(const cv::Mat &seq, const int max_iter, cv::Mat &TRANS, cv::Mat &EMIS, cv::Mat &INIT,bool UseUniformPrior = false)
{
    /* A Revealing Introduction to Hidden Markov Models, Mark Stamp */
    // 1. Initialization
    int iters = 0;
    int T = seq.cols; // number of element per sequence
    int C = seq.rows; // number of sequences
    int N = TRANS.rows; // number of states | also N = TRANS.cols | TRANS = A = {aij} - NxN
    int M = EMIS.cols; // number of observations | EMIS = B = {bj(k)} - NxM
    correctModel(TRANS,EMIS,INIT);
    cv::Mat FTRANS,FINIT,FEMIS;
    if (UseUniformPrior)
        getUniformModel(N,M,FTRANS,FEMIS,FINIT);
    else
    {
        FTRANS = TRANS.clone();
        FEMIS = EMIS.clone();
        FINIT = INIT.clone();
    }
    double logProb = -DBL_MAX;
    double oldLogProb;
    int data = 0;
    do {
        oldLogProb = logProb;

        // compute a0
        cv::Mat a(N,T,CV_64F);
        cv::Mat c(1,T,CV_64F); c.at<double>(0,0) = 0;
        for (int i=0;i<N;i++)
        {
            a.at<double>(i,0) = INIT.at<double>(0,i)*EMIS.at<double>(i,seq.at<int>(0,0));
            c.at<double>(0,0) += a.at<double>(i,0);
        }
        // scale the a0(i)
        c.at<double>(0,0) = 1/c.at<double>(0,0);
        for (int i=0;i<N;i++)
            a.at<double>(i,0) *= c.at<double>(0,0);

        // 2. The a-pass
        // compute at(i)
        for (int t=1;t<T;t++)
        {
            c.at<double>(0,t) = 0;
            for (int i=0;i<N;i++)
            {
                a.at<double>(i,t) = 0;
                for (int j=0;j<N;j++)
                    a.at<double>(i,t) += a.at<double>(i,t-1)*TRANS.at<double>(j,i);
                a.at<double>(i,t) = a.at<double>(i,t) * EMIS.at<double>(i,seq.at<int>(data,t));
                c.at<double>(0,t)+=a.at<double>(i,t);
            }
            // scale at(i)
            c.at<double>(0,t) = 1/c.at<double>(0,t);
            for (int i=0;i<N;i++)
                a.at<double>(i,t)=c.at<double>(0,t)*a.at<double>(i,t);
        }
        // 3. The B-pass
        cv::Mat b(N,T,CV_64F);
        // Let Bt-1(i) = 1 scaled by Ct-1
        for (int i=0;i<N;i++)
            b.at<double>(i,T-1) = c.at<double>(0,T-1);
        // B-pass
        for (int t=T-2;t>-1;t--)
            for (int i=0;i<N;i++)
            {
                b.at<double>(i,t) = 0;
                for (int j=0;j<N;j++)
                    b.at<double>(i,t) += TRANS.at<double>(i,j)*EMIS.at<double>(j,seq.at<int>(data,t+1))*b.at<double>(j,t+1);
                // scale Bt(i) with same scale factor as at(i)
                b.at<double>(i,t) *= c.at<double>(0,t);
            }
        // 4. Compute  Yt(i,j) and Yt(i)
        double denom;
        int index;
        cv::Mat YN(N,T,CV_64F);
        cv::Mat YNN(N*N,T,CV_64F);
        for (int t=0;t<T-1;t++)
        {
            denom = 0;
            for (int i=0;i<N;i++)
                for (int j=0;j<N;j++)
                    denom += a.at<double>(i,t)*TRANS.at<double>(i,j)*EMIS.at<double>(j,seq.at<int>(data,t+1))*b.at<double>(j,t+1);
            index = 0;
            for (int i=0;i<N;i++)
            {
                YN.at<double>(i,t) = 0;
                for (int j=0;j<N;j++)
                {
                    YNN.at<double>(index,t) = (a.at<double>(i,t)*TRANS.at<double>(i,j)*EMIS.at<double>(j,seq.at<int>(data,t+1))*b.at<double>(j,t+1))/denom;
                    YN.at<double>(i,t)+=YNN.at<double>(index,t);
                    index++;
                }
            }
        }
        // 5. Re-estimate A,B and pi
        // re-estimate pi
        for (int i=0;i<N;i++)
            INIT.at<double>(0,i) = YN.at<double>(i,0);
        // re-estimate A
        double numer;
        index = 0;
        for (int i=0;i<N;i++)
            for (int j=0;j<N;j++)
            {
                numer = 0;
                denom = 0;
                for (int t=0;t<T-1;t++)
                {
                    numer += YNN.at<double>(index,t);
                    denom += YN.at<double>(i,t);
                }
                TRANS.at<double>(i,j) = numer/denom;
                index++;
            }
        // re-estimate B
        for (int i=0;i<N;i++)
            for (int j=0;j<M;j++)
            {
                numer = 0;
                denom = 0;
                for (int t=0;t<T-1;t++)
                {
                    if (seq.at<int>(data,t)==j)
                        numer+=YN.at<double>(i,t);
                    denom += YN.at<double>(i,t);
                }
                EMIS.at<double>(i,j) = numer/denom;
            }
        correctModel(TRANS,EMIS,INIT);
        FTRANS = (FTRANS*(data+1)+TRANS)/(data+2);
        FEMIS = (FEMIS*(data+1)+EMIS)/(data+2);
        FINIT = (FINIT*(data+1)+INIT)/(data+2);
        // 6. Compute log[P(O|y)]
        logProb = 0;
        for (int i=0;i<T;i++)
            logProb += log(c.at<double>(0,i));
        logProb *= -1;
        // 7. To iterate or not
        data++;
        if (data >= C)
        {
            data = 0;
            iters++;
        }
    } while (iters<max_iter && logProb>oldLogProb);
    correctModel(FTRANS,FEMIS,FINIT);
    TRANS = FTRANS.clone();
    EMIS = FEMIS.clone();
    INIT = FINIT.clone();
}





int main()
{
    std::cout << "First we define Transition, Emission\nand Initial Probabilities of the model\n\n";
    double TRANSdata[] = {0.5 , 0.5 , 0.0,
                          0.0 , 0.7 , 0.3,
                          0.0 , 0.0 , 1.0};
    cv::Mat TRANS = cv::Mat(3,3,CV_64F,TRANSdata).clone();
    double EMISdata[] = {2.0/4.0 , 2.0/4.0 , 0.0/4.0 , 0.0/4.0 ,
                         0.0/4.0 , 2.0/4.0 , 2.0/4.0 , 0.0/4.0 ,
                         0.0/4.0 , 0.0/4.0 , 2.0/4.0 , 2.0/4.0 };
    cv::Mat EMIS = cv::Mat(3,4,CV_64F,EMISdata).clone();
    double INITdata[] = {1.0  , 0.0 , 0.0};
    cv::Mat INIT = cv::Mat(1,3,CV_64F,INITdata).clone();
//    CvHMM hmm;
//    hmm.printModel(TRANS,EMIS,INIT);
    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------
    std::cout << "\nAs an example, we generate 25 sequences each with 20 observations\nper sequence using the defined Markov model\n";
    srand ((unsigned int) time(NULL) );
    cv::Mat seq,states;
    generate(20,25,TRANS,EMIS,INIT,seq,states);
    std::cout << "\nGenerated Sequences:\n";
    for (int i=0;i<seq.rows;i++)
    {
        std::cout << i << ": ";
        for (int j=0;j<seq.cols;j++)
            std::cout << seq.at<int>(i,j);
        std::cout << "\n";
    }
    std::cout << "\nGenerated States:\n";
    for (int i=0;i<seq.rows;i++)
    {
        std::cout << i << ": ";
        for (int j=0;j<seq.cols;j++)
            std::cout << states.at<int>(i,j);
        std::cout << "\n";
    }
    std::cout << "\n";
    //----------------------------------------------------------------------------------

    std::cout << "\nProblem 3: Given an observation sequence O (can be several observations),\n";
    std::cout << "how do we find a model that maximizes the probability of O ?\n";
    std::cout << "The answer is by using the Baum-Welch algorithm to train a model.\n";
    std::cout << "To demonstrate this, initially we define a model by guess\n";
    std::cout << "and we estimate the parameters of the model for all the sequences\n";
    std::cout << "that we already got.\n";
    double TRGUESSdata[] = {2.0/3.0 , 1.0/3.0 , 0.0/3.0,
                            0.0/3.0 , 2.0/3.0 , 1.0/3.0,
                            0.0/3.0 , 0.0/3.0 , 3.0/3.0};
    cv::Mat TRGUESS = cv::Mat(3,3,CV_64F,TRGUESSdata).clone();
    double EMITGUESSdata[] = {1.0/4.0 , 1.0/4.0 , 1.0/4.0 , 1.0/4.0 ,
                              1.0/4.0 , 1.0/4.0 , 1.0/4.0 , 1.0/4.0 ,
                              1.0/4.0 , 1.0/4.0 , 1.0/4.0 , 1.0/4.0 };
    cv::Mat EMITGUESS = cv::Mat(3,4,CV_64F,EMITGUESSdata).clone();
    double INITGUESSdata[] = {0.6  , 0.2 , 0.2};
    cv::Mat INITGUESS = cv::Mat(1,3,CV_64F,INITGUESSdata).clone();

    train(seq,100,TRGUESS,EMITGUESS,INITGUESS);
    printModel(TRGUESS,EMITGUESS,INITGUESS);
    //----------------------------------------------------------------------------------
    std::cout << "\ndone.\n";
    return 0;
}
