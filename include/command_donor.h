# ifndef __COMMAND_DONOR__H
# define __COMMAND_DONOR_H
# include <command.h>

class CommandDonor :public Command
{
	private:
		Team	*donorTeam;
		Team	*acceptorTeam;
		void 	donorFinder();
		void	acceptorFinder();
		int	donor_or_acceptor_first;
	public:
		CommandDonor();
		void	setFlag(int Flag);
		void	setTeams(Team *a,Team *b);
		virtual void Run();
};
# endif
