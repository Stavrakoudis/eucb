# ifndef __COMMAND_RDF__H
# define __COMMAND_RDF__H
# include <command.h>

class CommandRdf: public Command
{
	private:
		int rdf_id;
		string rdf_criteria;
		Data radious;
		Data radious_step;
		Data radious_max;
		Data  volume;
		string	xst_filename;
		Data 	xst_freq;
	public:
		CommandRdf(int id,string criteria);
		void	setRadious(Data r1,Data r2,Data r3);
		void	setVolume(Data v);
		void	setXstFilename(string s);
		void	setXstFreq(Data x);
		virtual void Run();
		~CommandRdf();
};

# endif
