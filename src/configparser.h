/**
 *
 *  @file configparser.h
 *  C++ Declaration: ConfigParser
 *  @author Lorenz Quack <l_quac01@uni-muenster.de>
 *  Created on: 2009
 *  @author Georg Bergner <g.bergner@uni-muenster.de>
 *  2011: removed some bugs and unwanted features.
 *  Copyright: See COPYING file that comes with this distribution
 *  Taken from the supersymmetric Yang-Mills simulation program.
 *  Functions from the cpp-file made inline to avoid additional compilation units.
 */

#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

#include <exception>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <iomanip>
#include <iostream>
#include <algorithm>

namespace Parser {
typedef std::map<const std::string, std::string> section_t; ///< Collection of options in a section.
typedef std::map<const std::string, section_t> sections_t; ///< In this form the file is saved in the config parser.

/**
 General exception used by the ConfigParser class.
 */

class BaseError: public std::exception {

public:
	BaseError(const std::string &msg);
	virtual ~BaseError() throw () {
	}

	virtual const char *what() const throw ();

protected:
	std::string message;
};

/**
 Exception used by the ConfigParser class.
 */

class NoSectionError: public BaseError {

public:
	NoSectionError(const std::string &section);
	virtual ~NoSectionError() throw () {
	}

	const std::string section;
};

/**
 Exception used by the ConfigParser class.
 */

class DuplicateSectionError: public BaseError {

public:
	DuplicateSectionError(const std::string &section);
	virtual ~DuplicateSectionError() throw () {
	}

	const std::string section;
};

/**
 Exception used by the ConfigParser class.
 */

class NoOptionError: public BaseError {

public:
	NoOptionError(const std::string &option, const std::string &section);
	virtual ~NoOptionError() throw () {
	}

	const std::string section;
	const std::string option;
};

/**
 Exception used by the ConfigParser class.
 */

class ParsingError: public BaseError {

public:
	ParsingError(const std::string &filename);
	virtual ~ParsingError() throw () {
	}

	void append(const int lineno, const std::string &line);
	const std::string filename;
	std::vector<std::pair<int, std::string> > errors;
};

/**
 Exception used by the ConfigParser class.
 */

class MissingSectionHeaderError: public BaseError {

public:
	MissingSectionHeaderError(const std::string &filename, const int lineno,
			const std::string &line);
	virtual ~MissingSectionHeaderError() throw () {
	}

	const std::string filename;
	const int lineno;
	const std::string line;
};

/**
 Exception used by the ConfigParser class.
 */

class TypeConversionError: public BaseError {

public:
	TypeConversionError(const std::string &section, const std::string &option,
			const std::string &value, const std::string &type);
	virtual ~TypeConversionError() throw () {
	}

	const std::string section;
	const std::string option;
	const std::string value;
	const std::string type;
};

/**
 *  This class reads config files similar to the config parser  of the Python Standard Library.
 *  The ConfigParser class implements a basic configuration file parser
 *  language which provides a structure similar to what you would find on
 *  Microsoft Windows INI files.
 *  The configuration file consists of sections, led by a "[section]" header
 *  and followed by "name: value" entries, "name=value" is also accepted.
 *  Note that leading whitespace is removed from values.
 *  Lines beginning with '#' or ';' are ignored and may be used to provide
 *  comments.
 *
 *  This is more or less a C++ port of the RawConfigParser from the
 *  ConfigParser module from the Python Standard Library. However, not all
 *  features of that module have been implemented, yet.
 *  written by Lorenz Quack (l_quac01@uni-muenster.de)
 *  Reference: http://docs.python.org/library/configparser.html
 *
 *  Deleted some features of the basic implementation by Lorenz Quack and resolved some bugs
 *  Georg Bergner <G.Bergner@uni-muenster.de>
 *  @author Lorenz Quack <l_quac01@uni-muenster.de>
 *                                                                           */

class ConfigParser {

public:
	ConfigParser() throw ();
	ConfigParser(const ConfigParser& configParser) throw ();
	virtual ~ConfigParser() throw ();
	ConfigParser& operator=(const ConfigParser &configParser) throw (); // not sure if copying realy needed
	/**
	 * parse a config file filename.
	 * @param filename parsed file
	 * @return true for success false otherwise.
	 */
	bool read(const std::string& filename);

	/**
	 * parse a list of config files. (not realy tested -- not really needed)
	 * @param filenames string of filenames
	 * @return vector of files successfully read.
	 */
	std::vector<std::string> read(const std::vector<std::string>& filenames);
	/**
	 * Read a specific option in a section.
	 * @param section section of the value. [section]
	 * @param option option of the value. option = value
	 * @return value as a string.
	 */
	std::string get(const std::string& section, const std::string& option) const;
	/**
	 *  Like #get() but converts string to int.
	 * Throws excetion if conversion error.
	 */
	int getint(const std::string& section, const std::string& option) const;
	/**
	 *  Like #get() but converts string to a std vector of int until the first not interpretable character.
	 *  Throws excetion if vector empty.
	 */
	std::vector<int> getintvector(const std::string& section,
			const std::string& option) const ;
	/**
	 *  Like #get() but converts string to float.
	 * Throws excetion if conversion error.
	 */
	float getfloat(const std::string& section, const std::string& option) const;
	/**
	 *  Like #get() but converts string to double.
	 * Throws excetion if conversion error.
	 */
	double getdouble(const std::string& section,
			const std::string& option) const;
	/**
	 *  Like #get() but converts string to a std vector of double until the first not interpretable character.
	 *  Throws excetion if vector empty.
	 */
	std::vector<double> getdoublevector(const std::string& section,
			const std::string& option) const;
	/**
	 *  Like #get() but converts string to bool.
	 * accepts 0 or 1 or true or false as value
	 * Throws excetion if conversion error.
	 */
	bool getboolean(const std::string& section, const std::string& option) const;
	/**
	 * Returns true if section is present.
	 */
	bool has_section(const std::string& section) const throw ();
	/**
	 * Returns true if option in section is present.
	 */
	bool has_option(const std::string& section,
			const std::string& option) const throw ();
	/**
	 * Returns raw data: of all sections with their options.
	 */
	const sections_t& sections() const throw ();
	/**
	 * Returns list of all sections.
	 */
	std::vector<std::string> section_names() const throw ();
	/** returns a std::vector with all option names in section. */
	std::vector<std::string> option_names(const std::string& section) const;
	/** returns a map\<std::string, std::string\> which maps option names
	 to the option values of section. */
	const section_t& items(const std::string& section) const;
	/** sets the \<option\> in section to value.
	 If section doesn't exist a NoSectionError is thrown. */
	template<typename T>
	void set(const std::string& section, const std::string& option,
			const T& value);
	/** adds section to the configuration. Note that you need to
	 call "write" to save changes to disk. */
	void add_section(const std::string& section);
	/** attempts to remove section. returns "true" if section
	 existed "false" otherwise. */
	bool remove_section(const std::string& section) throw ();
	/** attempts to remove option from section.
	 throws a "NoSectionError" if section doesn't exist.
	 returns "true" if option existed "false" otherwise. */
	bool remove_option(const std::string& section, const std::string& option);
	/** writes the configuration to a file called filename */
	void write(const std::string& filename) const;

	/** overwrite this to control case sensitivity.
	 the default is case-insesitive.
	 */
	std::string optionxform(const std::string &optionstr) const throw ();

	/**
	 * Template version of #get() with conversion to type T.
	 * @param section section of value.
	 * @param option option of value.
	 * @return value converted to T.
	 */
	template<class T>
	T get_value(const std::string& section, const std::string& option) const{
		T ret;
		const std::string opt = this->get(section, option);

		if (from_string<T>(&ret, opt))
			return (ret);
		else
			throw TypeConversionError(section, option, opt, "template");
	}
  template<class T>
  T get_value_noerr(const std::string& section, const std::string& option, const T& defaultValue) const{
    T ret(defaultValue);
    try{
      ret=get_value<T>(section,option);
    }catch(...){
      return defaultValue;
    }
    return ret;
  }

private:
	/** the actuall workhorse of the parsing process. */
	void read_(std::ifstream& f, const std::string filename);
	/** useful function to convert a string str to a given type T and
	 * store it in t.
	 * returns true if the conversion was successful and false otherwise.
	 */
	template<typename T>
	static bool from_string(T *t, const std::string &str,
			std::ios_base & (*f)(std::ios_base&) =std::dec) {
		std::istringstream iss(str);
		return (!(iss >> f >> *t).fail());
	}

	/** writes the configuration to a std::string. */
	std::string write_() const throw ();
	sections_t *sections_; ///< raw data.
};

template<>
inline bool ConfigParser::from_string<std::vector<int> >(std::vector<int>* t,
		const std::string &str, std::ios_base & (*f)(std::ios_base&)) {
	std::istringstream iss(str);
	int num(0);
	for (iss >> num; iss.good(); iss >> num)
		t->push_back(num);
	if (iss.eof())
		t->push_back(num); // no conversion error but end of string is reached.
	return (!t->empty());
}

template<>
inline bool ConfigParser::from_string<std::vector<unsigned int> >(
		std::vector<unsigned int>* t, const std::string &str,
		std::ios_base & (*f)(std::ios_base&)) {
	std::istringstream iss(str);
	unsigned int num(0);
	for (iss >> num; iss.good(); iss >> num)
		t->push_back(num);
	if (iss.eof())
		t->push_back(num); // no conversion error but end of string is reached.
	return (!t->empty());
}

template<>
inline bool ConfigParser::from_string<std::vector<double> >(
		std::vector<double>* t, const std::string &str,
		std::ios_base & (*f)(std::ios_base&)) {
	std::istringstream iss(str);
	double num(0);
	for (iss >> num; iss.good(); iss >> num)
		t->push_back(num);
	if (iss.eof())
		t->push_back(num); // no conversion error but end of string is reached.
	return (!t->empty());
}

// template -> explicit inline.
template<typename T>
inline void ConfigParser::set(const std::string &section,
		const std::string &option, const T &value){
	sections_t::iterator sectionsIt = this->sections_->find(section);

	if (sectionsIt == sections_->end())
		throw NoSectionError(section);

	std::ostringstream tmp;

	tmp.setf(std::ios::boolalpha);

	tmp << value;

	(sectionsIt->second)[option] = tmp.str();
}

/* Here all inlining*/
inline ConfigParser::ConfigParser() throw () :
			sections_(NULL) {
	this->sections_ = new sections_t;
}

inline ConfigParser::ConfigParser(const ConfigParser& configParser) throw () :
			sections_(NULL) {
	this->sections_ = new sections_t;
	(*this->sections_) = configParser.sections();
}

inline ConfigParser::~ConfigParser() throw () {
	delete sections_;
}

inline ConfigParser& ConfigParser::operator=(
		const ConfigParser &configParser) throw () {
	delete sections_;
	this->sections_ = new sections_t;
	(*this->sections_) = configParser.sections();
	return (*this);
}

inline bool ConfigParser::read(const std::string& filename){
	std::ifstream f(filename.c_str());

	if (f.good()) {
		try {
			this->read_(f, filename);
		} catch (...) {
			std::cerr
					<< "\nERROR: An unhandled error occurred during parsing of config file \"%s\".\n"
					<< std::endl;
			f.close();
			throw;
		}

		f.close();

		return (true);
	}

	return (false);
}

inline std::vector<std::string> ConfigParser::read(
		const std::vector<std::string>& filenames) {
	std::vector<std::string> read_ok;

	for (std::vector<std::string>::const_iterator it = filenames.begin();
			it != filenames.end(); ++it) {
		if (read(*it))
			read_ok.push_back(*it);
	}

	return (read_ok);
}

inline const sections_t& ConfigParser::sections() const throw () {
	return (*this->sections_);
}

inline std::string ConfigParser::get(const std::string& section,
		const std::string& option) const {
	sections_t::iterator sectionIt = this->sections_->find(section);
	section_t::iterator optionIt;

	if (sectionIt == this->sections_->end())
		throw NoSectionError(section);
	else if ((optionIt = sectionIt->second.find(option))
			!= sectionIt->second.end())
		return (optionIt->second);
	else
		throw NoOptionError(option, section);
}

inline int ConfigParser::getint(const std::string& section,
		const std::string& option) const {
	int ret;
	const std::string opt = this->get(section, option);

	if (from_string<int>(&ret, opt))
		return (ret);
	else
		throw TypeConversionError(section, option, opt, "int");
}

inline std::vector<int> ConfigParser::getintvector(const std::string& section,
		const std::string& option) const  {
	std::vector<int> ret;
	const std::string opt = this->get(section, option);

	if (from_string<std::vector<int> >(&ret, opt))
		return (ret);
	else
		throw TypeConversionError(section, option, opt, "int");
}

inline float ConfigParser::getfloat(const std::string &section,
		const std::string &option) const  {
	float ret;
	const std::string opt = this->get(section, option);

	if (from_string<float>(&ret, opt))
		return (ret);
	else
		throw TypeConversionError(section, option, opt, "float");
}

inline double ConfigParser::getdouble(const std::string &section,
		const std::string &option) const {
	double ret;
	const std::string opt = this->get(section, option);

	if (from_string<double>(&ret, opt))
		return (ret);
	else
		throw TypeConversionError(section, option, opt, "double");
}

inline std::vector<double> ConfigParser::getdoublevector(
		const std::string& section, const std::string& option) const {
	std::vector<double> ret;
	const std::string opt = this->get(section, option);

	if (from_string<std::vector<double> >(&ret, opt))
		return (ret);
	else
		throw TypeConversionError(section, option, opt, "double");
}

inline bool ConfigParser::getboolean(const std::string &section,
		const std::string &option) const{
	bool ret;
	const std::string opt = this->get(section, option);

	if (from_string<bool>(&ret, opt, std::boolalpha))
		return (ret);
	else if (from_string<bool>(&ret, opt))
		return (ret);
	else
		throw TypeConversionError(section, option, opt, "bool");
}

inline void ConfigParser::add_section(const std::string &section){
	if (this->sections_->find(section) != this->sections_->end())
		throw DuplicateSectionError(section);

	(*this->sections_)[section];
}

inline bool ConfigParser::remove_section(const std::string &section) throw () {
	sections_t::iterator it = this->sections_->find(section);

	if (it == this->sections_->end())
		return (false);

	this->sections_->erase(it);

	return (true);
}

inline bool ConfigParser::remove_option(const std::string &section,
		const std::string &option)  {
	sections_t::iterator sectionsIt = this->sections_->find(section);

	if (sectionsIt == sections_->end())
		throw NoSectionError(section);

	section_t::iterator sectionIt = sectionsIt->second.find(option);

	if (sectionIt == sectionsIt->second.end())
		return (false);

	sectionsIt->second.erase(sectionIt);

	return (true);
}

inline void ConfigParser::write(const std::string &filename) const{
	std::ofstream f(filename.c_str());

	if (f.fail()) {
		f.close();
		std::stringstream msg;
		msg << "could not open " << filename << " for writing.";
		throw std::ios_base::failure(msg.str());
	}

	f << this->write_();

	f.close();
}

inline std::string ConfigParser::write_() const throw () {
	std::ostringstream buffer;

	for (sections_t::const_iterator it2 = sections_->begin();
			it2 != sections_->end(); ++it2) {
		buffer << "[" << it2->first << "]" << std::endl;

		for (section_t::const_iterator it = it2->second.begin();
				it != it2->second.end(); ++it) {
			// for each line break insert a tab character to make line continuations work
			std::string value = it->second;

			for (size_t pos = value.find("\n"); pos != std::string::npos; pos =
					value.find("\n", pos + 1))
				value.replace(pos, 1, "\n\t");

			buffer << it->first << " = " << value << std::endl;
		}

		buffer << std::endl;
	}

	return (buffer.str());
}

inline void ConfigParser::read_(std::ifstream &f, const std::string filename) {
	section_t* currentSection = NULL;
	int lineno = 0;
	std::string whitespaces = " \t\n\r\f\v";
	std::string optname, value;
	bool isExceptionPending = false;
	ParsingError pendingException(filename);

	while (true) {
		size_t pos1, pos2;
		char buffer[1024];

		if (!f.getline(buffer, 1024))
			break;

		std::string line = buffer;

		lineno++;

		/* ignore empty lines and comments */
		if ((line.find_first_not_of(whitespaces) == std::string::npos)
				|| (line.compare(0, 1, "#") == 0))
			continue;

		/* ignore lines starting with "rem" (BASIC-style comment) */
		if ((line.find("rem") == 0) || (line.find("REM") == 0)
				|| (line.find("Rem") == 0))
			continue;

		pos1 = line.find_last_not_of(whitespaces);

		line = line.substr(0, pos1 + 1); // strip of last whitepsaces.

		pos1 = line.find("[");

		pos2 = line.rfind("]");

		/* is it a section header? */
		if ((pos1 != std::string::npos) && (pos2 != std::string::npos)
				&& (pos1 < pos2)) {
			std::string sectionName = line.substr(pos1 + 1, pos2 - pos1 - 1);
			currentSection = &(*this->sections_)[sectionName]; // This either appends a secion or returns pointer to an existing section.
			optname.clear();
		} else if (currentSection == NULL) /* no section header in the file? */
			throw(MissingSectionHeaderError(filename, lineno, line));
		else /* an option line? */
		{
			pos1 = line.find_first_of(":=");

			if (pos1 != std::string::npos) {
				pos2 = line.find_first_not_of(whitespaces + ":=", pos1); // begins search from pos1 forward
				/* strip whitespaces of option name */
				pos1 = line.find_last_not_of(whitespaces, pos1 - 1); //begins search form pos1-1 backward
				optname = line.substr(0, pos1 + 1);
				/* check if there is a value */

				if (pos2 != std::string::npos) {
					/* strip whitespaces of option value */
					value = line.substr(pos2, std::string::npos);
					pos2 = value.find(";");

					if (pos2 != std::string::npos) {
						/* ";" is a comment delimiter only if it follows
						 *  a spacing character
						 */
						pos1 = value.find_first_not_of(whitespaces);
						pos2 = value.find_last_not_of(whitespaces, pos2 - 1);
						value = value.substr(pos1, pos2 - pos1 + 1);
					}

					/* allow empty values */
					if (value.compare("\"\"") == 0)
						value = "";
				} else {
					value = "";
				}

				(*currentSection)[optname] = value;
			} else {
				/* a non-fatal parsing error occurred.  set up the
				 * exception but keep going. the exception will be
				 * thrown at the end of the file and will contain a
				 * list of all bogus lines
				 */
				if (!isExceptionPending)
					isExceptionPending = true;

				pendingException.append(lineno, line);
			}
		}
	}

	/* if any parsing errors occurred, raise an exception */
	if (isExceptionPending)
		throw(pendingException);
}

inline std::vector<std::string> ConfigParser::section_names() const throw () {
	std::vector<std::string> sectionList;

	for (sections_t::iterator it = this->sections_->begin();
			it != this->sections_->end(); ++it)
		sectionList.push_back(it->first);

	return (sectionList);
}

inline std::vector<std::string> ConfigParser::option_names(
		const std::string &section) const{
	if (!has_section(section))
		throw(NoSectionError(section));

	std::vector<std::string> optionList;

	section_t* _section = &(*this->sections_)[section];

	for (section_t::iterator it = _section->begin(); it != _section->end();
			++it)
		optionList.push_back(it->first);

	return (optionList);
}

inline const section_t& ConfigParser::items(const std::string &section) const{
	if (!has_section(section))
		throw(NoSectionError(section));

	return ((*this->sections_)[section]);
}

inline std::string ConfigParser::optionxform(
		const std::string &optionstr) const throw () {
	std::string transformedOptionstr;
	std::transform(optionstr.begin(), optionstr.end(),
			std::back_inserter(transformedOptionstr), (int (*)(int) ) tolower);

return(	transformedOptionstr);
}

/* Exceptions */

inline BaseError::BaseError(const std::string &msg) :
			message(msg) {
}

inline const char *BaseError::what() const throw () {
	return (message.c_str());
}

inline NoSectionError::NoSectionError(const std::string &section) :
			BaseError(std::string("No section: ") + section),
			section(section) {
}

inline DuplicateSectionError::DuplicateSectionError(const std::string &section) :
			BaseError(std::string("No section: ") + section),
			section(section) {
}

inline NoOptionError::NoOptionError(const std::string &option,
		const std::string &section) :
			BaseError(
					std::string("No option ") + option + " in section: "
							+ section),
			section(section),
			option(option) {
}

inline ParsingError::ParsingError(const std::string &filename) :
			BaseError(std::string("File cpmtains parsing errors: ") + filename),
			filename(filename),
			errors() {
}

inline void ParsingError::append(const int lineno, const std::string &line) {
	std::stringstream msg;
	msg << std::endl << "\t[line " << std::setw(2) << lineno << "]: " << line;
	this->message += msg.str();
	this->errors.push_back(std::pair<int, std::string>(lineno, line));
}

inline MissingSectionHeaderError::MissingSectionHeaderError(
		const std::string &filename, const int lineno, const std::string &line) :
			BaseError(""),
			filename(filename),
			lineno(lineno),
			line(line) {
	std::stringstream msg;
	msg << "File contains no section headers." << std::endl << "file: "
			<< filename << ", line: " << lineno << std::endl << line;
	this->message = msg.str();
}

inline TypeConversionError::TypeConversionError(const std::string &section,
		const std::string &option, const std::string &value,
		const std::string &type) :
			BaseError(""),
			section(section),
			option(option),
			value(value),
			type(type) {
	std::stringstream msg;
	msg << "Could not convert option " << option << ": '" << value
			<< "' in section " << section << " to type: " << type << std::endl;
	this->message = msg.str();
}

inline bool ConfigParser::has_section(const std::string &section) const throw () {
	return (bool(this->sections_->find(section) != this->sections_->end()));
}

inline bool ConfigParser::has_option(const std::string &section,
		const std::string &option) const throw () {
	if (this->sections_->find(section) == this->sections_->end())
		return (false);
	else
		return (bool(
				(*sections_)[section].find(option)
						!= (*sections_)[section].end()));
}
} // namespace Parser
#endif /* CONFIGPARSER_H */
