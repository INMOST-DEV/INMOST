#ifndef INMOST_XML_INCLUDED
#define INMOST_XML_INCLUDED

#include "inmost.h"

namespace INMOST
{  
  std::string CharToHex(char c);
#if defined(USE_MESH)
  std::string ReferenceToString(INMOST::HandleType h, int pos);
#endif
#if defined(USE_AUTODIFF)
  std::string VariableToString(INMOST::Storage::var v);
#endif
  char * sstrip(char * str);
  std::string sstrip(const std::string & input);
  int ConvertHex(char in);
  char atoc(const char * str);
    

  class XMLReader
  {
    class Interpreter
    {
      bool error_state;
      std::vector<std::string> Expand(const std::string & input) const;
      std::vector<std::string> MakePolish(const std::vector<std::string> & input);
      void Print(const std::vector<std::string> & polish) const;
      double Run(const std::vector<std::string> & polish);
    public:
      Interpreter();
      Interpreter(const Interpreter & b);
      Interpreter & operator = (Interpreter const & b);
      double Evaluate(const std::string & str);
      bool isError();
      void ClearError();
    };
    Interpreter intrp;
    struct Stream
    {
      std::string src;
      std::istream * s;
      int linebreak, linechar;
      int hadlinebreak, hadlinechar;
    };
    std::vector<Stream> inp;
    int verbose;
    enum State
    {
      Intro, //read tag or read tag contents
      WaitTag, //wait tag name or comment
      ReadTag, //reading in tag name
      ReadCommentExclamation, //skipping comments
      ReadCommentQuestion, //skipping comments
      WaitAttribute, //reading attribute name
      ReadAttribute,
      WaitAttributeValue, //read attribute value
      ReadAttributeValue,
      ReadAttributeValueQuote,
      EndTag, //tag ended read closing
      WaitCloseTag,
      ReadCloseTagSlash,
      ReadCloseTagName,

      WaitContentsOpen,
      WaitContents,
      ReadContents, //parse a word
      ReadContentsVector, // {0,1,2,3}
      ReadContentsQuotes, // "hello world"
      ReadContentsMultiplier, // 123*5
      ReadContentsMultiplierSkopes, // 123*(SetSize/2)
      EndContents,
      ReadVector,

      EndOfFile, // end of file reached
      Failure //unexpected error
    } _state;
    Stream & get_Stream();
    const Stream & get_Stream() const;
    std::istream & get_iStream();
    const std::istream & get_iStream() const;
    //should not share the reference to the stream with another reader
    XMLReader(const XMLReader & other);
    XMLReader & operator =(XMLReader & other);
    char GetChar();
    //return one character back to the stream
    void RetChar();
    void SkipComments(State RetState);
    std::string StateName(State s) const;
  public:
  
    void Report(const char * fmt, ...) const;
    XMLReader(std::string sourcename, std::istream & input);
    void PushStream(std::string file);
    void PopStream();
    //read in <TagName returns TagName
    std::string ReadOpenTag();
    //read > or /> skipping for attributes
    int ReadCloseTag();
    bool isTagFinish() const;
    //read </TagName> or fail
    bool ReadFinishTag(std::string TagName);
    //read next attribute name, check isTagEnded
    std::string AttributeName();
    //read value of the attribute after reading it's name
    std::string AttributeValue();
    // > or /> was reached, should close ReadCloseTag
    // to finalize
    bool isTagEnded() const;
    //read in <![CDATA[
    //Quick and dirty, rewrite with states!
    bool ReadOpenContents();
    //get next full word inside contents
    std::string GetContentsWord();
    //read ]]>
    bool ReadCloseContents();
    //]]> was reached, should call ReadCloseContents
    //to finalize
    bool isContentsEnded() const;
    bool isFailure() const;
    bool isEof() const;
    
#if defined(USE_MESH)
    INMOST::ElementType atoes(const char * _str);
    INMOST::ElementType atoe(const char * _str);
    std::pair<INMOST::ElementType,int> atoh(const char * _str);
    std::pair<std::string,std::pair<INMOST::ElementType,int> > atorh(const char * _str);
#endif
#if defined(USE_AUTODIFF)
    INMOST::Storage::var atov(const char * _str);
#endif
    int EvaluateExpression(std::string expression);
    int ConvertMultiplier(std::string expression, int SetSize);
    void SplitValueMultiplier(std::string expression, std::string & value, std::string & multiplier);
    bool ParseBool(std::string word);
    void ParseCommaSeparated(std::string word, std::vector<std::string> & parsed, char symbol = ',');
    void ParseReal(std::string word, std::vector<INMOST::Storage::real> & Vector, int & Repeat, int SetSize);
#if defined(USE_AUTODIFF)
    void ParseVariable(std::string word, std::vector<INMOST::Storage::var> & Vector, int & Repeat, int SetSize);
#endif
    void ParseInteger(std::string word, std::vector<INMOST::Storage::integer> & Vector, int & Repeat, int SetSize);
    void ParseBulk(std::string word, std::string & Vector, int & Repeat, int SetSize);
    void ParseReference(std::string word, std::vector<std::pair<INMOST::ElementType,int> > & Vector, int & Repeat, int SetSize);
    void ParseRemoteReference(std::string word, std::vector< std::pair<std::string,std::pair<INMOST::ElementType,int> > > & Vector, int & Repeat, int SetSize);
    
    struct XMLAttrib
    {
      std::string name;
      std::string value;
    };

    struct XMLTag
    {
      std::string name; //<Name of the XML tag
      std::vector<XMLAttrib> attributes; //<List of attributes
      int finish; //<Whether to close the tag
      ///Was not able to read the tag
      bool Failure() const;
      ///Tag was red, can process the contents
      bool Process() const;
      ///Tag was not red, finish of enclosing tag was encountered
      bool Finalize() const;
      const XMLAttrib & GetAttib(int n) const;
      XMLAttrib & GetAttib(int n);
      int NumAttrib() const;
    };

    XMLTag OpenTag();
    bool CloseTag(XMLTag & tag);
  };
}


#endif //INMOST_XML_INCLUDED
