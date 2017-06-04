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
  std::string VariableToString(INMOST::variable v);
#endif
  
    

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
	/// 2 - lot's of output
	/// 1 - output key steps, currently in ReadXML
	void SetVerbosity(int verbosity) {verbose =  verbosity;}
    void Report(const char * fmt, ...) const;
    XMLReader(std::string sourcename, std::istream & input);
    void PushStream(std::string file);
    void PopStream();
	//wait for '<' on input,
	//returns true and changes state to WaitTag if '<' encountered,
	//otherwise returns false
	bool ExpectOpenTag();
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
    INMOST::variable atov(const char * _str);
#endif
    int EvaluateExpression(std::string expression);
    int ConvertMultiplier(std::string expression, int SetSize);
    void SplitValueMultiplier(std::string expression, std::string & value, std::string & multiplier);
    bool ParseBool(std::string word);
    void ParseCommaSeparated(std::string word, std::vector<std::string> & parsed, char symbol = ',');
#if defined(USE_MESH)
    void ParseReal(std::string word, std::vector<INMOST::Storage::real> & Vector, int & Repeat, int SetSize);
#if defined(USE_AUTODIFF)
    void ParseVariable(std::string word, std::vector<INMOST::Storage::var> & Vector, int & Repeat, int SetSize);
#endif
    void ParseInteger(std::string word, std::vector<INMOST::Storage::integer> & Vector, int & Repeat, int SetSize);
#endif
    void ParseBulk(std::string word, std::string & Vector, int & Repeat, int SetSize);
#if defined(USE_MESH)
    void ParseReference(std::string word, std::vector<std::pair<INMOST::ElementType,int> > & Vector, int & Repeat, int SetSize);
    void ParseRemoteReference(std::string word, std::vector< std::pair<std::string,std::pair<INMOST::ElementType,int> > > & Vector, int & Repeat, int SetSize);
#endif
    
	/// Structure for xml attribute.
    struct XMLAttrib
    {
      std::string name; //< Name of the attribute.
      std::string value; //< Value of the attribute.
    };

	/// Structure for xml tag with attributes.
    struct XMLTag
    {
      std::string name; //<Name of the XML tag.
      std::vector<XMLAttrib> attributes; //<List of attributes.
      int finish; //<Whether to close the tag.


	  ///This is data without ![CDATA[ wrap.
	  bool RawData() const {return finish == 5;}
	  ///This is data within ![CDATA[ wrap.
	  bool BlockData() const {return finish == 4;}
      ///Was not able to read the tag.
      bool Failure() const {return finish == 0;}
      ///Tag was red and have internal contents, can process the contents.
      bool Process() const {return finish == 1;}
	  ///Tag was red but do not have internal contents.
	  bool Stub() const {return finish == 2;}
      ///Tag was not red, finish of enclosing tag was encountered.
      bool Finalize() const {return finish == 3;}
	  ///Retrive attribute number n.
      const XMLAttrib & GetAttib(int n) const {return attributes[n];}
	  ///Retrive attribute number n.
      XMLAttrib & GetAttrib(int n) {return attributes[n];}
	  ///Retrive number of attributes.
      int NumAttrib() const {return (int)attributes.size();}
	  ///Retrive the name of the tag.
	  std::string GetName() const {return name;}
    };

    XMLTag OpenTag();
    bool CloseTag(XMLTag & tag);

	/// Structure defining entire XML file.
	struct XMLTree
	{
		XMLTag tag; //< tag information, such as name and attributes.
		std::vector<XMLTree> children; //< Children inside XML tag.
		std::string contents; //< Text inside of the tag.

		///Return next occurance of XML tag with the specified
		///name. Returns NumChildren() if not found.
		int FindChild(std::string name, int offset = -1) const;
		///Return next occurance of XML attribute with the specified
		///name. Returns NumAttrib() if not found.
		int FindAttrib(std::string name, int offset = -1) const;
		///Retrive a child of current XML tag with number n.
		const XMLTree & GetChild(int n) const {return children[n];}
		///Retrive a child of current XML tag with name
		///Returns NULL if not found.
		const XMLTree * GetChild(std::string name) const;
		///Retrive a child of current XML tag with attribute
		///Returns NULL if not found.
		const XMLTree * GetChildWithAttrib(std::string name, std::string value) const;
		///Retrive number of children.
		int NumChildren() const {return (int)children.size();}
		///Retrive attribute of current XML tag with number n.
		const XMLAttrib & GetAttrib(int n) const {return tag.GetAttib(n);}
		///Retrive attribute of current XML tag with name.
		///Returns NULL if not found.
		const std::string & GetAttrib(std::string name) const;
		///Retrive number of attributes.
		int NumAttrib() const {return tag.NumAttrib();}
		///Retrive the name of the tag.
		std::string GetName() const {return tag.GetName();}
		///Retrive contents of the tag.
		const std::string & GetContents() const {return contents;}
	};


private:
	std::string ReadUntil(std::string stop);
	int ReadXMLSub(XMLTree & root);
public:
	/// Read entire XML file into structure,
	/// it may be more efficient to read the file incrementially, depending on the size.
	/// See mesh_xml_file.cpp for incremential read.
	XMLTree ReadXML();
  };

  void WriteXML(const XMLReader::XMLTree & t, std::ostream & output, int offset = 0);

	typedef std::vector<XMLReader::XMLTree>::iterator xml_reader_tree_iterator_t;
	typedef std::vector<XMLReader::XMLAttrib>::iterator xml_reader_attrib_iterator_t;
}


#endif //INMOST_XML_INCLUDED
