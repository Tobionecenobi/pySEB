#include "SymbolInterface.hpp"
#include <cctype>

#ifdef USE_GINAC_IMPL
#include "GiNaCExpression.hpp"
#endif

using namespace std;

/*Resets Symbol interface*/
SymbolInterface* SymbolInterface::myInstance = 0;

/*In first call to instance, generate symbol interface class*/
SymbolInterface* SymbolInterface::instance(){
    if( !myInstance ){
        myInstance = new SymbolInterface();
    }
    return myInstance;
}

/*
    Replaces all occurances with the from string to the to string in str.
    PROBLEM: we need case insensitive replace.
*/
void SymbolInterface::replaceAll(std::string& str, const std::string& from, const std::string& to) {
    if(from.empty()) return;
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}

/*Replace psi only is the first four characters are psi_*/
void SymbolInterface::replacePsi(std::string& str){
    //Check that "psi" is the first characters in the string and the next character is "_"
    if(str.substr(0,4) == "psi_" || str.substr(0,4) == "Psi_"){
        str.replace(0,4,"\\psi_");
    }
}

/*Replace beta only if the first five characters are beta_*/
void SymbolInterface::replaceBeta(std::string& str){
    //Check that "beta" is the first characters in the string and the next character is "_"
    if(str.substr(0,5) == "beta_" || str.substr(0,5) == "Beta_"){
        str.replace(0,5,"\\beta_");
    }
}

/*Replace Rg only if the first two characters are Rg*/
void SymbolInterface::replaceRg(std::string& str){
    //Check that "Rg" is the first characters in the string and the next character is "_"
    if(str.substr(0,3) == "Rg_" || str.substr(0,3) == "rg_"){
        str.replace(0,2,"R_g");
        str.insert(str.size(), "}");
    }
}

string SymbolInterface::string2latex(const string& s)
/*
   String is one of F_<string>, A_<string>,<refstring>, Psi_<string>,<refstring>,<refstring>
   Here we transform one of these strings into LaTeX formatted expressions.

   f_<string>   ->  F_{string}

*/
{
   string latex(s);

   //Find the first occurance of _
    size_t pos=latex.find("_");

   // Transforms greek symbols
   replacePsi(latex);
   replaceBeta(latex);
   replaceRg(latex);
   replaceAll(latex, "#", "\\#" );

   // replace subscript _<string> with _{string}
   size_t found=latex.find("_" );
   if (found!=std::string::npos)     // We have a subscript
       {
          replaceAll(latex, "_", "_{");
          latex+="}";
       }

   return latex;
}

/*   A-Za-z0-9 and # is allowed.*/
bool SymbolInterface::testnamestring(const string& s)
{
    size_t pos=0;

    while (pos<s.size()-1)   // Don't test trailing zero.
      {
        if(!isalnum(s[pos])) return false;
        pos++;
      }
   return true;
}

string SymbolInterface::invalid_name_message(const string& s)
{
    for (size_t pos = 0; pos + 1 < s.size(); ++pos)
    {
        unsigned char ch = static_cast<unsigned char>(s[pos]);
        if (!std::isalnum(ch)) {
            return "Invalid character '" + string(1, s[pos]) + "' at position " + std::to_string(pos) + 
                   ". Allowed characters are A-Z, a-z, and 0-9.";
        }
    }
    return "Name validation failed. Allowed characters are A-Z, a-z, and 0-9.";
}


/*Test strings which can contain A-Za-z0-9, :.# is also allowed */
bool SymbolInterface::teststring(const string& s)
{
    size_t pos=0;

    while (pos<s.size()-1)   // Don't test trailing zero.
      {
        if(!isalnum(s[pos]) && s[pos]!=':'  && s[pos]!='#' && s[pos]!='.') return false;
        pos++;
      }
   return true;
}

string SymbolInterface::invalid_reference_message(const string& s)
{
    for (size_t pos = 0; pos + 1 < s.size(); ++pos)
    {
        unsigned char ch = static_cast<unsigned char>(s[pos]);
        if (!std::isalnum(ch) && s[pos] != ':' && s[pos] != '#' && s[pos] != '.') {
            return "Invalid character '" + string(1, s[pos]) + "' at position " + std::to_string(pos) +
                   ". Allowed characters are A-Z, a-z, 0-9, ':', '.', and '#'.";
        }
    }
    return "Reference validation failed. Allowed characters are A-Z, a-z, 0-9, ':', '.', and '#'.";
}



/*
   Gets a symbol based on a supplied string. If the string is already known, return the
   symbol we created the last time.
*/

const SymExprPtr& SymbolInterface::get( const string& s ){

// Store lowercase name as key to retrieve symbol
    string ss(s);  //transform(ss.begin(), ss.end(), ss.begin(), ::tolower);

    auto it = symbolDirectory.find( ss );
    if( it != symbolDirectory.end()) return it -> second;
                                else {
        // Create a new symbol using the symbolic factory
        SymExprPtr sym = SymbolicExpression::symbol(s);
        return symbolDirectory.insert( make_pair( ss, sym )).first->second;
    }
}


// For instance Rg, L, X
Expression SymbolInterface::getSymbol( const string& s ){

    if (!testnamestring(s)) throw SEBException("Bad symbol in variable name:"+s,"getSymbol(const string&)");
    return Expression(get(s));
}

// For instance F  <subunit>
Expression SymbolInterface::getSymbol( const string& s, const string& tag ){

    if (!teststring(s))      throw SEBException("Bad symbol in variable name:"+s,"getSymbol(const string&,const string&)");
    if (!teststring(tag))    throw SEBException("Bad symbol in variable name:"+tag,"getSymbol(const string&,const string&)");

    string si = s+"_"+tag;
    return Expression(get(si));
}

// s_tag:ref
Expression SymbolInterface::getSymbol( const string& s, const string& tag, const string& ref ){

    if (!teststring(s))      throw SEBException("Bad symbol in variable name:"+s  , "getSymbol(const string&,const string&,const string&)");
    if (!teststring(tag))    throw SEBException("Bad symbol in variable name:"+s  , "getSymbol(const string&,const string&,const string&)");
    if (!teststring(ref))     throw SEBException("Bad symbol in index name:"+tag , "getSymbol(const string&,const string&,const string&)");

    string si;
    si = s+"_"+tag+":"+ref;
    return Expression(get(si));
}

//  E.g.  psi   sub_unit/structure   r1, r2
Expression SymbolInterface::getSymbol(const string& s, const string& tag, const string& index1, const string& index2 ){

    if (!teststring(s))      throw SEBException("Bad symbol in variable name:"+s,"getSymbol(const string&)");
    if (!teststring(tag))    throw SEBException("Bad symbol in variable name:"+tag  , "getSymbol(const string&,const string&,const string&)");
    if (!teststring(index1))  throw SEBException("Bad symbol in index name:"+index1 , "getSymbol(const string&,const string&,const string&)");
    if (!teststring(index2))  throw SEBException("Bad symbol in index name:"+index2 , "getSymbol(const string&,const string&,const string&)");

    string si;
    if (index1<index2) si = s+"_"+tag+":"+index1+","+index2;
                  else si = s+"_"+tag+":"+index2+","+index1;

    return Expression(get(si));
}

// Provide symbol table
map<string, SymExprPtr> SymbolInterface::getSymbolTable()
{
    return symbolDirectory;
}




string to_string_cform(const SymExprPtr& expr)
{
    return expr->to_cform();
}

string to_string_python(const SymExprPtr& expr)
{
    return expr->to_python();
}

string to_string_latex(const SymExprPtr& expr)
{
    return expr->to_latex();
}

string to_string(const SymExprPtr& expr)
{
    return expr->to_string();
}

// The implementation of SymbolicFactory methods is in SymbolicInterface.cpp



