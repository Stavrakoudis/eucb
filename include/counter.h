# ifndef __COUNTER__H
# define __COUNTER__H
# include <util.h>

typedef struct
{
	int		frame;
	int		total;
	vector<string> 	residue;
	IntVector	count;	
}COUNT;

class	Counter
{
	private:
		vector<COUNT>  count;
		vector<string> diff_res;
	public:
		Counter();
		int 	size();
		void	add(int Frame,string Residue,int Count);
		void	print(string filename);
		void	printStat(string filename);
		void	printHist(string filename);
		void	printSmooth(string filename);
};
# endif
