# ifndef __COMMANDENTROPY__H
# define __COMMANDENTROPY__H

# include <command.h>
# include <eigen.h>

class CommandEntropy : public Command
{
	private:
		int nframes;
		Data temp;
		int	perresidue;
		int	fitflag;
		int 	pairflag;
		void	getEntropies(int eigendim,DoubleVector &evalue,Data &entropy1,Data &entropy2);
		void	getSubEntropy(IntVector &posAtom,IntVector &starFrame,
					DoubleVector &entropy1,DoubleVector &entropy2);
	public:
		CommandEntropy(int n,Data t,int f,int fit,int pair);
		virtual void Run();
		~CommandEntropy();
};

# endif
