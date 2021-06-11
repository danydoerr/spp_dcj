/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

/*

This file contains XML util functions

Created the: 09-02-2015
by: Wandrille Duchemin

Last modified the: 02-09-2016
by: Wandrille Duchemin

*/

#include "XMLUtils.h"




/*
Takes:
 - line (string): a line out of a phyloxml file

Returns:
 (string): name of the balise in the name. "" if there is no balise 
*/
string getLineBaliseName(string line)
{
	size_t baliseBegin = line.find("<");

	bool Closing = false;

	string baliseName;

	if(baliseBegin == string::npos) // character not found -> ignore
		return "";

	baliseBegin++;
	
	if(line[baliseBegin] == '/')
	{
		Closing = true;
		baliseBegin++;
	}

	size_t EndBaliseName = string::npos;
	size_t spacePos = line.find(" ",baliseBegin);
	if(spacePos != string::npos)
		EndBaliseName = spacePos;
	size_t greaterPos = line.find(">",baliseBegin);
	if((greaterPos != string::npos) && (( greaterPos < EndBaliseName) || (EndBaliseName == string::npos)))
		EndBaliseName = greaterPos;
	if(EndBaliseName == string::npos)
		return "";
	else
		baliseName = line.substr(baliseBegin,EndBaliseName - baliseBegin);

	return baliseName;
}

/*
 - line (string ) : line in an XLM file
 - properties (map <string, string> * ) : pointer to a map where keys will be attribute name and associated values will be the value of these attribute
 - value (string * ) : eventual value between the tags

 Returns:
 	(string):

*/
string InterpretLineXML(string line, map <string, string> * properties, string * value)
{
	// 1. find the name of the balise

	size_t baliseBegin = line.find("<");

	string baliseName;

	if(baliseBegin == string::npos) // no balise
		return "";

	baliseBegin++;

	size_t EndBaliseName = string::npos;
	size_t spacePos = line.find(" ",baliseBegin);
	if(spacePos != string::npos)
		EndBaliseName = spacePos;
	size_t greaterPos = line.find(">",baliseBegin);
	if((greaterPos != string::npos) && (( greaterPos < EndBaliseName) || (EndBaliseName == string::npos)))
		EndBaliseName = greaterPos;


	if(EndBaliseName == string::npos)
		return "";// no balise
	else
		baliseName = line.substr(baliseBegin,EndBaliseName - baliseBegin);	

	//balise name found
	//2. look for properties
	size_t currentIndex = EndBaliseName;

	currentIndex++;

	while(currentIndex < greaterPos)
	{
		//reading a property
		string PropertyName = "";
		string PropertyValue = "";

		while(line[currentIndex] != '=')
		{
			PropertyName.push_back(line[currentIndex]);
			currentIndex++;
		}
		currentIndex++;
		if(line[currentIndex] == '"') // if first character of the value is a ", then we go to the next one
		{
			currentIndex++;
			while(line[currentIndex] != '"') // reading the value
			{
				PropertyValue.push_back(line[currentIndex]);
				currentIndex++;
			}
			currentIndex++;		
		}
		else // go to the next > or  ' '
		{
			while((line[currentIndex] != '>') && (line[currentIndex] != ' ') ) // reading the value
			{
				PropertyValue.push_back(line[currentIndex]);
				currentIndex++;
			}
		}
		//PropertyValue set, we are either on a ' ' or a '>'
		currentIndex++; // going over that character for next round

		(*properties)[PropertyName] = PropertyValue;

		//cout << "added property " <<  PropertyName << " -> " << PropertyValue << endl;
	}
	//we are after the > character

	size_t lowerpos = line.find("<",currentIndex);

	if(( lowerpos == string::npos)||( lowerpos == currentIndex)) //no value 
		return baliseName;

	//extracting value
	*value = *value + line.substr(currentIndex, lowerpos - currentIndex);

//	cout << line.substr(currentIndex, lowerpos - currentIndex) << endl;
//	cout << *value << endl;

	return baliseName;
}


/*
Takes:
	- ifstream& fileIN : file the data comes from
	- string tag : tag to find the end of
Returns:
	(bool) : true if the closing tag was found
*/
bool goToTagEnd(ifstream& fileIN , string tag)
{
	return goToNextTag(fileIN , "/" + tag );
}

/*
Takes:
	- ifstream& fileIN : file the data comes from
	- string tag : tag to find the next of

Returns:
	(bool) : true if the next tag was found
*/
bool goToNextTag(ifstream& fileIN , string tag)
{
	if( fileIN.eof() )
		return false;

	string line;
	getline( fileIN, line );

	map <string, string> * properties = new map <string,string>;
	string * value = new string(""); 

	string Tname = InterpretLineXML(line,properties,value);

	

	while(Tname.compare(tag) != 0) 
	{
		if(fileIN.eof())
			return false;

		getline( fileIN, line );
		Tname = InterpretLineXML(line,properties,value);
	}

	return true;
}