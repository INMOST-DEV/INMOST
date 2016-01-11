#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "inmost.h"
#include <stdarg.h> //for va_list

#if defined(USE_MESH)



std::string CharToHex(char c)
{
  char const hex[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B','C','D','E','F'};
  std::string str = "";
  str.append(&hex[(c & 0xF0) >> 4], 1);
  str.append(&hex[c & 0xF], 1);
  return str;
}

std::string ReferenceToString(INMOST::HandleType h, int pos)
{
  std::stringstream ret;
  switch(INMOST::GetHandleElementType(h))
  {
  case INMOST::NODE: ret << "Node:"; break;
  case INMOST::EDGE: ret << "Edge:"; break;
  case INMOST::FACE: ret << "Face:"; break;
  case INMOST::CELL: ret << "Cell:"; break;
  case INMOST::ESET: ret << "Set:"; break;
  case INMOST::MESH: ret << "Mesh:"; break;
  }
  ret << pos;//INMOST::GetHandleID(h);
  return ret.str();
}



#if defined(USE_AUTODIFF)
std::string VariableToString(INMOST::Storage::var v)
{
  std::stringstream ret;
  const INMOST::Sparse::Row & r = v.GetRow();
  ret << "(";
  ret << v.GetValue() << ";";
  ret << r.Size();
  if( !r.Empty() )
  {
    ret << ";";
    for(int q = 0; q < (int)r.Size()-1; ++q)
    {
      ret << r.GetValue(q) << ";";
      ret << r.GetIndex(q) << ";";
    }
    ret << r.GetValue(r.Size()-1) << ";";
    ret << r.GetIndex(r.Size()-1);
  }
  ret << ")";
  return ret.str();
}
#endif

class XMLReader
{
  class Interpreter
  {
    bool error_state;
    std::vector<std::string> Expand(std::string & input)
    {
      std::vector<std::string> ret;
      std::string put;
      for(int k = 0; k < (int)input.length(); ++k)
      {
        if( get_priority(input[k]) != -1 )
        {
          if( !put.empty() ) 
          {
            ret.push_back(put);
            put.clear();
          }
          ret.push_back(std::string(1,input[k]));
        }
        else put.push_back(input[k]);
      }
      if( !put.empty() ) ret.push_back(put);
      return ret;
    }
    int get_priority(char c)
    {
      switch(c)
      {
      case '(': return 0;
      case ')': return 1;
      case '+':
      case '-': return 8;
      case '*':
      case '/': return 9;
      case '~': return 10;
      default: return -1;
      }
    }
    std::vector<std::string> MakePolish(std::vector<std::string> & input)
    {
      std::vector<std::string> stack, ret;
      int priority;
      for(int k = 0; k < (int)input.size(); ++k)
      {
        if( input[k].size() == 1 && (priority = get_priority(input[k][0])) != -1 )
        {
          char op = input[k][0];
          if( op != '(' )
          {
            while(!stack.empty() && priority <= get_priority(stack.back()[0]) )
            {
              ret.push_back(stack.back());
              stack.pop_back();
            }
          }
          if( op == ')' )
          {
            if( stack.empty() )
            {
              std::cout << "Warning: parentheses unbalanced" << std::endl;
              error_state = true;
              break;
            }
            else if( stack.back()[0] == '(' )
              stack.pop_back();
            else 
            {
              std::cout << "Warning: expected left bracket" << std::endl;
              error_state = true;
              break;
            }
          }
          else if( op == '-' && (k == 0 || (input[k-1].size() == 1 && get_priority(input[k-1][0]) != -1)) ) //unary minus
            stack.push_back("~");
          else
            stack.push_back(input[k]);
        }
        else ret.push_back(input[k]);
      }
      while(!stack.empty())
      {
        ret.push_back(stack.back());
        stack.pop_back();
      }
      return ret;
    }
    void Print(std::vector<std::string> & polish)
    {
      for(int k = 0; k < (int)polish.size(); ++k)
        std::cout << polish[k] << " ";
      std::cout << std::endl;
    }
    double Run(std::vector<std::string> & polish)
    {
      std::vector<double> stack;
      for(int k = 0; k < (int)polish.size(); ++k)
      {
        if( polish[k].size() == 1 && get_priority(polish[k][0]) != -1 )
        {
          double larg, rarg;
          char op = polish[k][0];
          switch(op)
          {
          case '+':
            if( stack.size() < 2 )
            {
              std::cout << "Less then two arguments in stack for + operand" << std::endl;
              error_state = true;
              return 1.0e+20;
            }
            rarg = stack.back();
            stack.pop_back();
            larg = stack.back();
            stack.pop_back();
            stack.push_back(larg+rarg);
            break;
          case '-':
            if( stack.size() < 2 )
            {
              std::cout << "Less then two arguments in stack for - operand" << std::endl;
              error_state = true;
              return 1.0e+20;
            }
            rarg = stack.back();
            stack.pop_back();
            larg = stack.back();
            stack.pop_back();
            stack.push_back(larg-rarg);
            break;
          case '*':
            if( stack.size() < 2 ) 
            {
              std::cout << "Less then two arguments in stack for * operand" << std::endl;
              error_state = true;
              return 1.0e+20;
            }
            rarg = stack.back();
            stack.pop_back();
            larg = stack.back();
            stack.pop_back();
            stack.push_back(larg*rarg);
            break;
          case '/':
            if( stack.size() < 2 ) 
            {
              std::cout << "Less then two arguments in stack for / operand" << std::endl;
              error_state = true;
              return 1.0e+20;
            }
            rarg = stack.back();
            stack.pop_back();
            larg = stack.back();
            stack.pop_back();
            stack.push_back(larg/rarg);
            break;
          case '~':
            if( stack.size() < 1 ) 
            {
              std::cout << "No arguments in stack for unary minus operand" << std::endl;
              error_state = true;
              return 1.0e+20;
            }
            larg = stack.back();
            stack.pop_back();
            stack.push_back(-larg);
            break;
          }
        }
        else stack.push_back(atof(polish[k].c_str()));
      }
      if( stack.size() != 1 )
      {
        std::cout << "There are more operands on stack, but no operators" << std::endl;
        error_state = true;
        return 1.0e+20;
      }
      return stack.back();
    }
  public:
    Interpreter() :error_state(false) {}
    Interpreter(const Interpreter & b)
      : error_state(b.error_state)
    {}
    Interpreter & operator = (Interpreter const & b)
    {
      error_state = b.error_state;
      return * this;
    }
    double Evaluate(std::string & str)
    {
      //const char * debug_str = str.c_str();
      std::vector<std::string> decompose = Expand(str);
      std::vector<std::string> polish = MakePolish(decompose);
      //Print(polish);
      return Run(polish);
    }
    bool isError() {return error_state;}
    void ClearError() {error_state = false;}
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
  Stream & get_Stream() {return inp.back();}
  std::istream & get_iStream() {return *inp.back().s;}
  //should not share the reference to the stream with another reader
  XMLReader(const XMLReader & other) {}
  XMLReader & operator =(XMLReader & other) {return *this;}
  char GetChar()
  {
    char c = '\0';
    get_iStream().get(c);
    get_Stream().hadlinebreak = get_Stream().linebreak;
    get_Stream().hadlinechar = get_Stream().linechar;
    if( c == '\n' ) 
    {
      ++get_Stream().linebreak;
      get_Stream().linechar = 0;
    }
    else ++get_Stream().linechar;
    if( get_iStream().eof() ) 
    {
      if( inp.size() > 1 )
      {
        PopStream();
        c = GetChar();
      }
      else
        _state = EndOfFile;
    }
    if( get_iStream().fail() )
    {
      Report("Stream failed while getting the char");
      _state = Failure;
    }
    return c;
  }
  //return one character back to the stream
  void RetChar()
  {
    get_Stream().linebreak = get_Stream().hadlinebreak;
    get_Stream().linechar = get_Stream().hadlinechar;
    get_iStream().unget();
    if( get_iStream().fail() ) 
    {
      Report("Stream failed while ungetting the char");
      _state = Failure;
    }
  }

  void SkipComments(State RetState)
  {
    int ntmp = 0;
    char tmp[3] = {'\0','\0','\0'};
    char c;
    bool done = false;
    if( _state == ReadCommentExclamation )
    {
      while(!done)
      {
        c = GetChar();
        tmp[ntmp] = c;
        if( tmp[ntmp] == '>' && tmp[(ntmp-1+3)%3] == '-' && tmp[(ntmp-2+3)%3] == '-' ) 
        {
          _state = RetState;
          done = true;
        }
        ntmp = (ntmp+1)%3;
      }
    }
    else if( _state == ReadCommentQuestion )
    {
      while(!done)
      {
        c = GetChar();
        tmp[ntmp] = c;
        if( tmp[ntmp] == '>' && tmp[(ntmp-1+2)%2] == '?' ) 
        {
          _state = RetState;
          done = true;
        }
        ntmp = (ntmp+1)%2;
      }
    }
    else 
    {
      Report("Unexpected state %s while reading comments",StateName(_state).c_str());
      _state = Failure; //What are we doing here?
    }
  }
  std::string StateName(State s)
  {
    switch(s)
    {
    case Intro: return "Intro";
    case WaitTag: return "WaitTag";
    case ReadTag: return "ReadTag";
    case ReadCommentExclamation: return "ReadCommentExclamation";
    case ReadCommentQuestion: return "ReadCommentQuestion";
    case WaitAttribute: return "WaitAttribute";
    case ReadAttribute: return "ReadAttribute";
    case WaitAttributeValue: return "WaitAttributeValue";
    case ReadAttributeValue: return "ReadAttributeValue";
    case ReadAttributeValueQuote: return "ReadAttributeValueQuote";
    case EndTag: return "EndTag";

    case ReadVector: return "ReadVector";

    case WaitContentsOpen:  return "WaitContentsOpen";
    case WaitContents:  return "WaitContents";
    case ReadContents: return "ReadContents";
    case ReadContentsVector: return "ReadContentsVector";
    case ReadContentsMultiplier: return "ReadContentsMultiplier";
    case ReadContentsMultiplierSkopes: return "ReadContentsMultiplierSkopes";
    case ReadContentsQuotes: return "ReadContentsQuotes";
    case EndContents: return "EndContents";

    case WaitCloseTag: return "WaitCloseTag";
    case ReadCloseTagSlash: return "WaitCloseTagSlash";
    case ReadCloseTagName: return "WaitCloseTagName";
    case EndOfFile: return "EndOfFile";
    case Failure: return "Failure";
    };
    return "Unspecified";
  }
public:
  
  void Report(const char * fmt, ...) 
  { 
    std::cout << get_Stream().src << ":row:" << get_Stream().linebreak << ":col:" << get_Stream().linechar << " ";
    {
      char stext[16384];
	    va_list ap;
      if ( fmt == NULL ) {std::cout << std::endl; return;}
	    va_start(ap,fmt);
	    vsprintf(stext,fmt,ap);
	    va_end(ap);
      std::cout << stext;
    }
    std::cout << std::endl;
  }
  XMLReader(std::string sourcename, std::istream & input) :intrp(),_state(Intro)
  {
    Stream add;
    add.src = sourcename;
    add.linebreak = 0;
    add.linechar = 0;
    add.hadlinebreak = 0;
    add.hadlinechar = 0;
    add.s = &input;
    inp.push_back(add); verbose = 0;
    if( get_iStream().fail() )
    {
      Report("Got a bad stream on input in %s",__FUNCTION__);
    }
  }
  void PushStream(std::string file) 
  {
    Stream add;
    add.linebreak = 0;
    add.linechar = 0;
    add.hadlinebreak = 0;
    add.hadlinechar = 0;
    add.src = file;
    add.s = new std::fstream(file.c_str(),std::ios::in);
    inp.push_back(add);
    if( get_iStream().fail() )
    {
      Report("Got a bad stream on input in %s",__FUNCTION__);
    }
  }
  void PopStream() 
  {
    if( inp.size() > 1 )
      delete static_cast<std::fstream *>(inp.back().s);
    inp.pop_back();
    if( _state == EndOfFile && !inp.empty() )
      _state = Intro;
  }
  //read in <TagName returns TagName
  std::string ReadOpenTag()
  {
    std::string ret;
    char c;
    bool done = false;
    if( !(_state == Intro) )
    {
      Report("Cannot open tag from state %s",StateName(_state).c_str());
      _state = Failure;
      return "";
    }
    while(!done)
    {
      c = GetChar(); 
      switch(_state)
      {
      case Intro:
        if( c == '<' ) 
        {
          if( verbose ) Report("info: waiting tag name");
          _state = WaitTag;
        }
        else if( !isspace(c) ) //do not expect anything except for spacing
        {
          Report("Unexpected text character %c",c);
          _state = Failure;
          done = true;
        }
        break;
      case WaitTag:
        if( c == '?' ) 
        {
          if( verbose ) Report("info: skipping comments");
          _state = ReadCommentQuestion;
          SkipComments(WaitTag);
        }
        else if( c == '!' ) 
        {
          if( verbose ) Report("info: skipping comments");
          _state = ReadCommentExclamation;
          SkipComments(WaitTag);
        }
        else if( c == '/' )
        {
          RetChar();
          _state = ReadCloseTagSlash;
          done = true;
        }
        else if( isalpha(c) )
        {
          if( verbose ) Report("info: reading tag name");
          ret.push_back(c);
          _state = ReadTag;
        }
        break;
      case ReadTag:
        if( isspace(c) ) 
        {
          if( verbose ) Report("info: waiting attribute name");
          done = true;
          _state = WaitAttribute;
        }
        else if( c == '/' || c == '>' )
        {
          if( verbose ) Report("info: tag ended");
          RetChar(); //push character back to the stream
          done = true;
          _state = EndTag;
        }
        else if( isalpha(c) ) ret.push_back(c);
        else Report("unexpected character %c in XML tag name",c);
        break;
      case EndOfFile: Report("Unexpected end of file while reading XML tag name"); done = true; break;
      case Failure: Report("Unrecoverable error while reading XML tag name"); done = true; break;
      default: Report("Unexpected state %s",StateName(_state).c_str()); done = true; break;
      }
    }
    if( verbose ) Report("info: opened tag %s",ret.c_str());
    return ret;
  }
  //read > or /> skipping for attributes
  int ReadCloseTag()
  {
    char tmp[2];
    tmp[0] = GetChar();
    if( tmp[0] == '>' )
    {
      _state = Intro;
      if( verbose ) Report("info: closed tag");
      return 1; //tag was finished with >
    }
    else if( tmp[0] == '/' ) //close single stage tag
    {
      tmp[1] = GetChar();
      if( tmp[1] == '>' )
      {
        _state = Intro;
        if( verbose ) Report("info: closed tag");
        return 2; //tag was halted with />
      }
      Report("Encountered '%c%c' while expecting '/>' for tag closing",tmp[0],tmp[1]);
    }
    Report("Encountered '%c' while expecting '>' for tag closing",tmp[0]);
    _state = Failure;
    return 0;
  }
  bool isTagFinish() {return _state == ReadCloseTagSlash;}
  //read </TagName> or fail
  bool ReadFinishTag(std::string TagName)
  {
    std::string name;
    bool done = false;
    char c;
    if( !(_state == Intro || _state == ReadCloseTagSlash) )
    {
      Report("Cannot read finish tag from state %s",StateName(_state).c_str());
      _state = Failure;
      return false;
    }
    if( _state == Intro ) _state = WaitCloseTag;
    while(!done)
    {
      c = GetChar();  
      switch(_state)
      {
      case WaitCloseTag:
        if( isspace(c) ) continue;
        if( c != '<' )
        {
          Report("Expected '<' instead of %c",c);
          _state = Failure;
          return false;
        }
        else _state = ReadCloseTagSlash;
        break;
      case ReadCloseTagSlash:
        if( c == '?' )
        {
          if( verbose ) Report("info: skipping comments");
          _state = ReadCommentQuestion;
          SkipComments(WaitCloseTag);
        }
        else if( c == '!' )
        {
          if( verbose ) Report("info: skipping comments");
          _state = ReadCommentExclamation;
          SkipComments(WaitCloseTag);
        }
        else if( c != '/' )
        {
          Report("Expected '/' instead of %c",c);
          _state = Failure;
          return false;
        }
        else _state = ReadCloseTagName;
        break;
      case ReadCloseTagName:
        if( isalpha(c) ) name.push_back(c);
        else if( c == '>' ) done = true;
        else Report("Unexpected symbol %c in tag name",c);
        break;
      case EndOfFile:
        Report("Unexpected end of file while searching for </%s>",TagName.c_str());
        return false;
      case Failure:
        Report("Unexpected failure while searching for </%s>",TagName.c_str());
        return false;
      default: Report("Unexpected state %s",StateName(_state).c_str()); return false;
              
      }
    }
    if( verbose ) Report("info: finished tag %s",name.c_str());
    _state = Intro;
    return name == TagName;
  }
  //read next attribute name, check isTagEnded
  std::string AttributeName()
  {
    std::string ret;
    bool done = false;
    char c;
    if( _state == EndTag ) return "";
    if( _state != WaitAttribute )
      Report("Attribute was not expected, state %s",StateName(_state).c_str());
    while(!done)
    {
      c = GetChar();
      switch(_state)
      {
      case WaitAttribute:
        if( isspace(c) ) continue;
        else if( c == '<' )
        {
          c = GetChar();
          if( c == '?' )
          {
            if( verbose ) Report("info: skipping comments");
            _state = ReadCommentQuestion;
            SkipComments(WaitAttribute);
          }
          else if( c == '!' )
          {
            if( verbose ) Report("info: skipping comments");
            _state = ReadCommentExclamation;
            SkipComments(WaitAttribute);
          }
          else 
          {
            Report("Expected a comment, got '<%c'",c);
            _state = Failure;
            done = true;
          }
        }
        else if( isalpha(c) ) 
        {
          if( verbose ) Report("info: reading attribute name");
          ret.push_back(c);
          _state = ReadAttribute;
        }
        else if( c == '>' || c == '/' )
        {
          if( verbose ) Report("info: tag ended");
          RetChar();
          done = true;
          _state = EndTag;
        }
        break;
      case ReadAttribute:
        if( isalpha(c) ) ret.push_back(c);
        else if( c == '=' || c == ' ' ) 
        {
          if( c == '=' ) RetChar();
          _state = WaitAttributeValue;
          done = true;
        }
        else 
        {
          Report("Unexpected symbol %c while reading attribute name",c);
          _state = Failure;
          done = true;
        }
        break;
      case EndOfFile:
        Report("Unexpected end of file while reading attribute name");
        done = true;
        break;
      case Failure:
        Report("Unexpected failure while reading attribute name");
        done = true;
        break;
      default: Report("Unexpected state %s",StateName(_state).c_str()); done = true; break;
      }
    }
    if( verbose ) Report("info: attribute name %s",ret.c_str());
    return ret;
  }
  //read value of the attribute after reading it's name
  std::string AttributeValue()
  {
    std::string ret;
    bool done = false;
    char c;
    if( _state == EndTag ) return "";
    if( _state != WaitAttributeValue )
      Report("Attribute was not expected, state %s",StateName(_state).c_str());
    while(!done)
    {
      c = GetChar();
      switch(_state)
      {
      case WaitAttributeValue:
        if( isspace(c) ) continue;
        else if( c == '=' )
        {
          if( verbose ) Report("info: reading attribute value");
          _state = ReadAttributeValue;
        }
        else if( c == '>' || c == '/' )
          Report("Unexpected end of XML tag while searching for '='");
        else Report("Unexpected character %c while searching for '='",c);
        break;
      case ReadAttributeValue:
        if( c == '"' && ret.empty() ) 
        {
          if( verbose ) Report("info: reading attribute value in quotes");
          _state = ReadAttributeValueQuote;
        }
        else if( c == '>' || c =='/' )
        {
          if( verbose ) Report("info: end of tag");
          _state = EndTag;
          done = true;
        }
        else if( !isspace(c) ) ret.push_back(c);
        else if( isspace(c) )
        {
          if( verbose ) Report("info: end reading attribute value");
          _state = WaitAttribute;
          done = true;
        }
        else Report("Unexpected symbol %c while reading attribute value",c);
        break;
      case ReadAttributeValueQuote:
        if( c == '"' )
        {
          if( verbose ) Report("info: end reading attribute value");
          _state = WaitAttribute;
          done = true;
        }
        else if( !isprint(c) )
        {
          Report("Unprintable character %x encountered",c);
          _state = Failure;
          done = true;
        }
        else ret.push_back(c);
        break;
      case EndOfFile:
        Report("Unexpected end of file while reading attribute name");
        done = true;
        break;
      case Failure:
        Report("Unexpected failure while reading attribute name");
        done = true;
        break;
      default: Report("Unexpected state %s",StateName(_state).c_str()); done = true; break;
      }
    }
    if( verbose ) Report("info: attribute value %s",ret.c_str());
    return ret;
  }
  // > or /> was reached, should close ReadCloseTag
  // to finalize
  bool isTagEnded()
  {
    return _state == EndTag;
  }
  //read in <![CDATA[
  //Quick and dirty, rewrite with states!
  bool ReadOpenContents()
  {
    std::string tmp;
    bool done = false;
    char c;
    if( _state != Intro )
    {
      Report("Cannot read contents opening from state %s",StateName(_state).c_str());
      _state = Failure;
      return false;
    }
    _state = WaitContentsOpen;
    while(!done)
    {
      c = GetChar();
      if( isspace(c) ) continue;
      else if( c == '<' || c == '!' || c == '[' || c == 'C' || c == 'D' || c == 'A' || c == 'T' || c == '?' || c == '-' )
        tmp.push_back(c);
      else Report("Unexpected character %c while reading '<![CDATA['",c);
      if( tmp.size() == 2 && tmp == "<?" )
      {
        tmp.clear();
        if( verbose ) Report("info: skipping comments");
        _state = ReadCommentQuestion;
        SkipComments(WaitContentsOpen);
      }
      else if( tmp.size() == 4 && tmp == "<!--" )
      {
        tmp.clear();
        if( verbose ) Report("info: skipping comments");
        _state = ReadCommentExclamation;
        SkipComments(WaitContentsOpen);
      }
      else if( tmp.size() == 9 )
      {
        if( tmp == "<![CDATA[" ) 
        {
          if( verbose ) Report("info: contents intro was read");
          _state = WaitContents;
          done = true;
        }
        else
        {
          Report("Expected '<![CDATA[' but got '%s'",tmp.c_str());
          _state = Failure;
          done = true;
        }
      }
    }
    return _state == WaitContents;
  }
  //get next full word inside contents
  std::string GetContentsWord()
  {
    std::string ret;
    char c;
    bool done = false;
    int testend = 0;
    if( _state == EndContents ) return "";
    if( _state != WaitContents )
    {
       Report("Cannot read contents from state %s",StateName(_state).c_str());
      _state = Failure;
      return "";
    }
    while(!done)
    {
      c = GetChar();
      /*
      if( !isprint(c) )
      {
        Report("Unprintable character %x on input",c);
        _state = Failure;
        break;
      }
      */
      if( c == ']' && !(_state == ReadContentsQuotes) )
        testend++;
      else if( c == '>' && testend == 2 )
      {
        if( ret[ret.size()-2] == ']' && ret[ret.size()-1] == ']' )
        {
          ret.resize(ret.size()-2);
          RetChar(); //return '>'
          _state = EndContents;
          break;
        }
        else testend = 0;
      }
      switch(_state)
      {
      case WaitContents:
        if( isspace(c) ) continue;
        if( c == '{' )
        {
          ret.push_back(c);
          _state = ReadContentsVector;
        }
        else if( c == '*' )
        {
          ret.push_back(c);
          _state = ReadContentsMultiplier;
        }
        else if( c == '"' )
        {
          ret.push_back(c);
          _state = ReadContentsQuotes;
        }
        else
        {
          ret.push_back(c);
          _state = ReadContents;
        }
        break;
      case ReadContents:
        if( isspace(c) )
        {
          done = true;
          _state = WaitContents;
        }
        else if( c == '*' )
        {
          ret.push_back(c);
          _state = ReadContentsMultiplier;
        }
        else ret.push_back(c);
        break;
      case ReadContentsVector:
        if( c == '}' )
        {
          ret.push_back(c);
          _state = ReadContents;
        }
        else if( !isspace(c) ) ret.push_back(c);
        break;
      case ReadContentsMultiplier:
        if( isspace(c) )
        {
          if( ret[ret.size()-1] != '*' ) //maybe user have put a space after the multiplier
          {
            done = true;
            _state = WaitContents;
          }
        }
        else if( c == '(' && ret[ret.size()-1] == '*' ) //expression for the skope
        {
          ret.push_back(c);
          _state = ReadContentsMultiplierSkopes;
        }
        else ret.push_back(c);
        break;
      case ReadContentsMultiplierSkopes:
        if( c == ')' ) //do not expect anything beyond that
        {
          ret.push_back(c);
          _state = WaitContents;
          done = true;
        }
        else if( !isspace(c) ) ret.push_back(c);
        break;
      case ReadContentsQuotes:
        if( c == '"' )
        {
          if( ret[ret.size()-1] == '\\' ) //this was a guarded quote
          {
            ret.resize(ret.size()-1); //remove guard
            ret.push_back(c); //put skope
          }
          else //may have a multiplier after the closing quote
          {
            ret.push_back(c);
            _state = ReadContents;
          }
        }
        else ret.push_back(c);
        break;
      case EndOfFile:
        Report("Unexpected end of file while reading attribute name");
        done = true;
        break;
      case Failure:
        Report("Unexpected failure while reading attribute name");
        done = true;
        break;
      default: Report("Unexpected state %s",StateName(_state).c_str()); done = true; break;
      }
    }
    return ret;
  }
  //read ]]>
  bool ReadCloseContents()
  {
    if( _state == EndContents && GetChar() == '>' )
    {
      _state = Intro;
      return true;
    }
    return false;
  }
  //]]> was reached, should call ReadCloseContents
  //to finalize
  bool isContentsEnded()
  {
    return _state == EndContents;
  }
  bool isFailure() {return _state == Failure;}
  bool isEof() {return _state == EndOfFile;}
  int ConvertHex(char in)
  {
    int ret = tolower(in) - 48;
    if( ret > 10 ) ret -= 7;
    if( ret > 15 ) ret -= 32;
    return ret;
  }
  char atoc(const char * str)
  {
    return (char)(ConvertHex(str[0])*16 + ConvertHex(str[1]));
  }
  INMOST::ElementType atoes(const char * _str)
  {
    INMOST::ElementType type = INMOST::NONE;
    std::string stype(_str);
    if( stype == "Cells" ) type = INMOST::CELL;
    else if( stype == "Edges" ) type = INMOST::EDGE;
    else if( stype == "Faces" ) type = INMOST::FACE;
    else if( stype == "Nodes" ) type = INMOST::NODE;
    else if( stype == "Sets" ) type = INMOST::ESET;
    else if( stype == "Mesh" ) type = INMOST::MESH;
    else if( stype == "None" ) type = INMOST::NONE;
    else Report("Unexpected keyword %s expected Nodes,Edges,Faces,Cells or Sets",_str);
    return type;
  }
  INMOST::ElementType atoe(const char * _str)
  {
    INMOST::ElementType type = INMOST::NONE;
    std::string stype(_str);
    if( stype == "Cell" ) type = INMOST::CELL;
    else if( stype == "Edge" ) type = INMOST::EDGE;
    else if( stype == "Face" ) type = INMOST::FACE;
    else if( stype == "Node" ) type = INMOST::NODE;
    else if( stype == "Set" ) type = INMOST::ESET;
    else if( stype == "Mesh" ) type = INMOST::MESH;
    else if( stype == "None" ) type = INMOST::NONE;
    else Report("Unexpected keyword %s expected Node,Edge,Face,Cell,Set or Mesh",_str);
    return type;
  }
  std::pair<INMOST::ElementType,int> atoh(const char * _str)
  {
    std::string str(_str);
    if( str == "None" ) return std::make_pair(INMOST::NONE,0);
    if( str == "None:0" ) return std::make_pair(INMOST::NONE,0);
    size_t dots = str.find(':');
    if( dots == std::string::npos ) return std::make_pair(INMOST::NONE,1);
    std::string stype = str.substr(0,dots);
    std::string offset = str.substr(dots+1,std::string::npos);
    INMOST::ElementType type = atoe(stype.c_str());
    int ioffset = atoi(offset.c_str());
    if( type == INMOST::NONE ) Report("Cannot understand element type %s",stype.c_str());
    return std::make_pair(type,ioffset);
  }
  std::pair<std::string,std::pair<INMOST::ElementType,int> > atorh(const char * _str)
  {
    std::string str(_str);
    if( str == "None" ) return std::make_pair("",std::make_pair(INMOST::NONE,0));
    if( str == ":None:0" ) return std::make_pair("",std::make_pair(INMOST::NONE,0));
    size_t dots1 = str.find(':');
    if( dots1 == std::string::npos ) return std::make_pair("",std::make_pair(INMOST::NONE,1));
    size_t dots2 = str.find(':',dots1+1);
    std::string meshname = str.substr(0,dots1);
    std::string stype = str.substr(dots1+1,dots2-dots1);
    std::string offset = str.substr(dots2+1,std::string::npos);
    INMOST::ElementType type = atoe(stype.c_str());
    int ioffset = atoi(offset.c_str());
    if( type == INMOST::NONE ) Report("Cannot understand element type %s",stype.c_str());
    return std::make_pair(meshname,std::make_pair(type,ioffset));
  }
#if defined(USE_AUTODIFF)
  INMOST::Storage::var atov(const char * _str)
  {
    INMOST::Sparse::Row entries;
    INMOST::Storage::real val;
    if( !(_str[0] == '(' && _str[1] == ')' ) ) Report("Expected scopes for variable %s",_str);
    std::string str(_str+1,strlen(_str)-2);
    std::vector<std::string> decomposed;
    ParseCommaSeparated(str,decomposed,';');
    if( decomposed.size() < 2 || decomposed.size()%2 != 0 ) Report("Malformed variable %s",_str);
    val = atof(decomposed[0].c_str());
    entries.Resize(atoi(decomposed[1].c_str()));
    if( decomposed.size()/2-1 != entries.Size() ) Report("Number of derivatives do not correspond to the number of entries %s",_str);
    for(int q = 2; q < (int)entries.Size(); q+=2)
    {
      entries.GetValue(q) = atof(decomposed[q].c_str());
      entries.GetIndex(q) = atoi(decomposed[q+1].c_str());
    }
    return INMOST::variable(val,entries);
  }
#endif
  int EvaluateExpression(std::string expression)
  {
    //port short RPN version later from DiscretizationToolkit/GPRS/interpreter.cpp
    return (int)intrp.Evaluate(expression);
    return atoi(expression.c_str()); //stupid version for now
  }
  char * sstrip(char * str)
  {
	  int text_start = 0, text_end = (int)strlen(str);
	  for(text_start = 0; isspace(str[text_start]) && text_start < text_end; text_start++);
	  if( text_start == text_end ) return str+text_start;
	  for(text_end = text_end-1; isspace(str[text_end]) && text_end > text_start; text_end--);
	  str[text_end+1] = '\0';
	  return str+text_start;
  }

  std::string sstrip(const std::string & input)
  {

	  char temp[2048];
	  strcpy(temp,input.c_str());
	  return std::string(sstrip(temp));
  }

  int ConvertMultiplier(std::string expression, int SetSize)
  {
    if( !expression.empty() )
    {
      //replace SetSize with value
      std::stringstream conv;
      conv << SetSize;
      size_t pos = 0;
      while((pos = expression.find("SetSize",pos)) != std::string::npos )
      {
        expression.replace(pos,7,conv.str());
        pos += conv.str().size();
      }
      return EvaluateExpression(expression);
    }
    else return 1;
  }
  void SplitValueMultiplier(std::string expression, std::string & value, std::string & multiplier)
  {
    size_t mult = expression.find('*');
    if( mult != std::string::npos )
    {
      value = expression.substr(0,mult);
      multiplier = expression.substr(mult+1,std::string::npos);
    }
    else
    {
      value = expression;
      multiplier = "";
    }
  }
  bool ParseBool(std::string word)
  {
    bool ret = false;
    if( word == "False" || word == "false" || word == "0" || word == "no" || word == "N" || word == "No" || word == "NO" ) ret = false;
    else if( word == "True" || word == "true" || word == "1" || word == "yes" || word == "Y" || word == "Yes" || word == "YES" ) ret = true;
    else Report("Cannot understand boolean value %s",word.c_str());
    return ret;
  }
  void ParseCommaSeparated(std::string word, std::vector<std::string> & parsed, char symbol = ',')
  {
    parsed.clear();
    size_t comma_prev = 0, comma;
    std::string substr;
    do
    {
      comma = word.find(symbol,comma_prev);
      substr = word.substr(comma_prev,comma-comma_prev);
      parsed.push_back(sstrip(substr));
      comma_prev = comma+1;
    } while( comma != std::string::npos );
  }
  void ParseReal(std::string word, std::vector<INMOST::Storage::real> & Vector, int & Repeat, int SetSize)
  {
    //const char * debug = word.c_str();
    std::string value, multiplier;
    SplitValueMultiplier(word,value,multiplier);
    Vector.clear();
    //const char * debug_value = value.c_str();
    //const char * debug_multiplier = multiplier.c_str();
    if( value[0] == '{' ) //multiple elements
    {
      size_t comma_prev = 1, comma;
      std::string substr;
      do
      {
        comma = value.find(',',comma_prev);
        if( comma == std::string::npos ) comma = value.find('}',comma_prev);
        substr = value.substr(comma_prev,comma-comma_prev);
        //const char * debug_substr = substr.c_str();
        Vector.push_back(atof(substr.c_str()));
        comma_prev = comma+1;
      } while( value[comma] != '}' );
    }
    else Vector.push_back(atof(value.c_str())); //single element

    Repeat = ConvertMultiplier(multiplier,SetSize);    
  }
#if defined(USE_AUTODIFF)
  void ParseVariable(std::string word, std::vector<INMOST::Storage::var> & Vector, int & Repeat, int SetSize)
  {
    std::string value, multiplier;
    SplitValueMultiplier(word,value,multiplier);
    Vector.clear();
    if( value[0] == '{' ) //multiple elements
    {
      size_t comma_prev = 1, comma;
      std::string substr;
      do
      {
        comma = value.find(',',comma_prev);
        if( comma == std::string::npos ) comma = value.find('}',comma_prev);
        substr = value.substr(comma_prev,comma-comma_prev);
        Vector.push_back(atov(substr.c_str()));
        comma_prev = comma+1;
      } while( value[comma] != '}' );
    }
    else Vector.push_back(atov(value.c_str())); //single element

    Repeat = ConvertMultiplier(multiplier,SetSize);    
  }
#endif
  void ParseInteger(std::string word, std::vector<INMOST::Storage::integer> & Vector, int & Repeat, int SetSize)
  {
    std::string value, multiplier;
    SplitValueMultiplier(word,value,multiplier);
    Vector.clear();
    if( value[0] == '{' ) //multiple elements
    {
      size_t comma_prev = 1, comma;
      std::string substr;
      do
      {
        comma = value.find(',',comma_prev);
        if( comma == std::string::npos ) comma = value.find('}',comma_prev);
        substr = value.substr(comma_prev,comma-comma_prev);
        Vector.push_back(atoi(substr.c_str()));
        comma_prev = comma+1;
      } while( value[comma] != '}' );
    }
    else Vector.push_back(atoi(value.c_str())); //single element

    Repeat = ConvertMultiplier(multiplier,SetSize);    
  }
  void ParseBulk(std::string word, std::string & Vector, int & Repeat, int SetSize)
  {
    std::string value, multiplier;
    SplitValueMultiplier(word,value,multiplier);
    Vector.clear();
    if( value[0] == '"' ) //just text
    {
      Vector.insert(0, value.c_str()+1, value.size()-2);
    }
    else if( value[0] == '{' ) //multiple elements
    {
      size_t comma_prev = 1, comma;
      std::string substr;
      do
      {
        comma = value.find(',',comma_prev);
        if( comma == std::string::npos ) comma = value.find('}',comma_prev);
        substr = value.substr(comma_prev,comma-comma_prev);
        if( substr.size() != 2 )
          Report("Unexpected size %d of substring %s in vector, expected 2",substr.size(),substr.c_str());
        Vector.push_back(atoc(substr.c_str()));
        comma_prev = comma+1;
      } while( value[comma] != '}' );
    }
    else if( value.size() >= 2 && value.size()%2 == 0 )
    {
      for(int k = 0; k < (int)value.size(); k+=2 )
        Vector.push_back(atoc(value.c_str()+k)); //single element
    }
    else Report("Unexpected size of input string %s size %d",value.c_str(),value.size());
    Repeat = ConvertMultiplier(multiplier,SetSize);
  }
  void ParseReference(std::string word, std::vector<std::pair<INMOST::ElementType,int> > & Vector, int & Repeat, int SetSize)
  {
    std::string value, multiplier;
    SplitValueMultiplier(word,value,multiplier);
    Vector.clear();
    if( value[0] == '{' )
    {
      size_t comma_prev = 1, comma;
      std::string substr;
      do
      {
        comma = value.find(',',comma_prev);
        if( comma == std::string::npos ) comma = value.find('}',comma_prev);
        substr = value.substr(comma_prev,comma-comma_prev);
        Vector.push_back(atoh(substr.c_str()));
        if( Vector.back().first == INMOST::NONE && Vector.back().second == 1 )
          Report("Cannot convert handle to the element, %s",substr.c_str());
        comma_prev = comma+1;
      } while( value[comma] != '}' );
    }
    else 
    {
      Vector.push_back(atoh(value.c_str())); 
      if( Vector.back().first == INMOST::NONE && Vector.back().second == 1 )
        Report("Cannot convert handle to the element, %s",value.c_str());
    }

    Repeat = ConvertMultiplier(multiplier,SetSize);
  }

  void ParseRemoteReference(std::string word, std::vector< std::pair<std::string,std::pair<INMOST::ElementType,int> > > & Vector, int & Repeat, int SetSize)
  {
    std::string value, multiplier;
    SplitValueMultiplier(word,value,multiplier);
    Vector.clear();
    if( value[0] == '{' )
    {
      size_t comma_prev = 1, comma;
      std::string substr;
      do
      {
        comma = value.find(',',comma_prev);
        if( comma == std::string::npos ) comma = value.find('}',comma_prev);
        substr = value.substr(comma_prev,comma-comma_prev);
        Vector.push_back(atorh(substr.c_str()));
        if( Vector.back().first == "" )
          Report("Cannot extract mesh name, %s",substr.c_str());
        if( Vector.back().second.first == INMOST::NONE && Vector.back().second.second == 1 )
          Report("Cannot convert handle to the element, %s",substr.c_str());
        comma_prev = comma+1;
      } while( value[comma] != '}' );
    }
    else 
    {
      Vector.push_back(atorh(value.c_str())); 
      if( Vector.back().first == "" )
          Report("Cannot extract mesh name, %s",value.c_str());
      if( Vector.back().second.first == INMOST::NONE && Vector.back().second.second == 1 )
        Report("Cannot convert handle to the element, %s",value.c_str());
    }

    Repeat = ConvertMultiplier(multiplier,SetSize);
  }

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
    bool Failure() {return finish == 0;}
    ///Tag was red, can process the contents
    bool Process() {return finish == 1 || finish == 2;}
    ///Tag was not red, finish of enclosing tag was encountered
    bool Finalize() {return finish == 3;}
    XMLAttrib & GetAttib(int n) {return attributes[n];}
    int NumAttrib() {return (int)attributes.size();}
  };

  XMLTag OpenTag()
  {
    std::string include = "";
    XMLTag ret;
    XMLAttrib attr;
    ret.name = ReadOpenTag();
    if( !isTagFinish() )
    {
      attr.name = AttributeName();
      while(!isTagEnded())
      {
        attr.value = AttributeValue();
        if( attr.name == "Include" )
          include = attr.value;
        else ret.attributes.push_back(attr);
        attr.name = AttributeName();
      }
      ret.finish = ReadCloseTag();
      if( !include.empty() ) PushStream(include.c_str());
    }
    else
    {
      ret.name = "";
      ret.finish = 3;
    }
    return ret;
  }
  bool CloseTag(XMLTag & tag)
  {
    if( tag.finish == 3 ) 
    {
      Report("%s:%d Trying to close the tag that is beyond the last tag",__FILE__,__LINE__);
      return false;
    }
    else if( tag.finish == 2 ) return true; //no need to read anything
    else if( tag.finish == 0 ) return false; //there was a Failure
    else return ReadFinishTag(tag.name);
  }
};

namespace INMOST
{
  void Mesh::LoadXML(std::string File)
  {
    std::fstream infile(File.c_str(),std::ios::in);
    XMLReader reader(File,infile);
    std::vector<Tag> tags;
    std::vector<HandleType> new_nodes, new_edges, new_faces, new_cells, new_sets;
    std::vector<ElementSet::ComparatorType> set_comparators;
    XMLReader::XMLTag PassTag;
    bool pass_tag = false;
    int nmeshes = 1;

    XMLReader::XMLTag TagParallelMesh = reader.OpenTag();
    if( TagParallelMesh.name != "ParallelMesh" )
    {
      reader.Report("Incorrect XML tag %s expected ParallelMesh",TagParallelMesh.name.c_str());
      throw BadFile;
    }

    for(int q = 0; q < TagParallelMesh.NumAttrib(); ++q)
    {
      XMLReader::XMLAttrib & attr = TagParallelMesh.GetAttib(q);
      if( attr.name == "Number" ) nmeshes = atoi(attr.value.c_str());
      else reader.Report("Unused attribute for ParallelMesh %s='%s'",attr.name.c_str(),attr.value.c_str());
    }
    
    for(XMLReader::XMLTag TagMesh = reader.OpenTag(); !TagMesh.Finalize() && TagMesh.name == "Mesh"; reader.CloseTag(TagMesh), TagMesh = reader.OpenTag())
    {
      bool repair_orientation = false;
      for(int q = 0; q < TagMesh.NumAttrib(); ++q)
      {
        XMLReader::XMLAttrib & attr = TagMesh.GetAttib(q);
        if( attr.name == "RepairOrientation" ) repair_orientation = reader.ParseBool(attr.value);
        else reader.Report("Unused attribute for Tags %s='%s'",attr.name.c_str(),attr.value.c_str());
      }

      { //Tags

        XMLReader::XMLTag TagTags = reader.OpenTag();
        
        if( TagTags.name != "Tags" )
        {
          reader.Report("Incorrect XML tag %s expected Tags",TagTags.name.c_str());
          throw BadFile;
        }

        int ntags = 0;
        bool matchntags = false;
        for(int q = 0; q < TagTags.NumAttrib(); ++q)
        {
          XMLReader::XMLAttrib & attr = TagTags.GetAttib(q);
          if( attr.name == "Number" ) 
          {
            ntags = atoi(attr.value.c_str());
            matchntags = true;
          }
          else reader.Report("Unused attribute for Tags %s='%s'",attr.name.c_str(),attr.value.c_str());
        }

        tags.reserve(ntags);


        XMLReader::XMLTag TagTag;
        for(TagTag = reader.OpenTag(); !TagTag.Finalize() && TagTag.name == "Tag"; reader.CloseTag(TagTag), TagTag = reader.OpenTag())
        {
          std::vector<std::string> parsed;
          std::string tagname = "";
          enumerator size = ENUMUNDEF;
          DataType type = DATA_REAL;
          ElementType sparse = NONE;
          ElementType defined = NONE;
          for(int q = 0; q < TagTag.NumAttrib(); ++q)
          {
            XMLReader::XMLAttrib & attr = TagTag.GetAttib(q);
            if( attr.name == "Name" ) tagname = attr.value;
            else if( attr.name == "Size" ) 
            {
              if( attr.value == "Variable" )
                size = ENUMUNDEF;
              else size = atoi(attr.value.c_str());
            }
            else if( attr.name == "Type" )
            {
              if( attr.value == "Real" ) type = DATA_REAL;
              else if( attr.value == "Integer" ) type = DATA_INTEGER;
              else if( attr.value == "Reference" ) type = DATA_REFERENCE;
              else if( attr.value == "RemoteReference" ) type = DATA_REMOTE_REFERENCE;
              else if( attr.value == "Bulk" ) type = DATA_BULK;
#if defined(USE_AUTODIFF)
              else if( attr.value == "Variable" ) type = DATA_VARIABLE;
#endif
            }
            else if( attr.name == "Sparse" )
            { 
              reader.ParseCommaSeparated(attr.value,parsed);
              for(int q = 0; q < (int)parsed.size(); ++q) sparse |= reader.atoes(parsed[q].c_str());
            }
            else if( attr.name == "Definition" )
            {
              reader.ParseCommaSeparated(attr.value,parsed);
              for(int q = 0; q < (int)parsed.size(); ++q) defined |= reader.atoes(parsed[q].c_str());
            }
            else reader.Report("Unused attribute for Tags %s='%s'",attr.name.c_str(),attr.value.c_str());
          }
          if( tagname == "" )
            reader.Report("Tag name was not specified");
          else if( defined == NONE )
            reader.Report("Domain of definition for the tag was not specified");
          tags.push_back(CreateTag(tagname,type,defined,sparse,size));
        }
        reader.CloseTag(TagTags);

        if( !TagTag.Finalize() )
        {
          PassTag = TagTag;
          pass_tag = true;
        }

        if( matchntags && ntags != tags.size() ) reader.Report("Number %d of XML tags Tag red do not match to the specified number %d",tags.size(),ntags);
      }

      { //Nodes
        int nnodes = 0, ndims = 3;
        XMLReader::XMLTag TagNodes;
        if( pass_tag )
        {
          TagNodes = PassTag;
          pass_tag = false;
        }
        else
          TagNodes = reader.OpenTag();

        if( TagNodes.name != "Nodes" )
        {
          reader.Report("Incorrect XML tag %s expected Nodes",TagNodes.name.c_str());
          throw BadFile;
        }

        for(int q = 0; q < TagNodes.NumAttrib(); ++q)
        {
          XMLReader::XMLAttrib & attr = TagNodes.GetAttib(q);
          if( attr.name == "Number" ) nnodes = atoi(attr.value.c_str());
          else if( attr.name == "Dimensions" )
          {
            ndims = atoi(attr.value.c_str());
            if( GetDimensions() != ndims ) SetDimensions(ndims);
          }
          else reader.Report("Unused attribute for Tags %s='%s'",attr.name.c_str(),attr.value.c_str());
        }

        new_nodes.reserve(nnodes);

        {
          if( reader.ReadOpenContents() )
          {
            std::vector<double> Vector;
            int Repeat;
            dynarray<Storage::real,3> xyz;
            for(std::string val = reader.GetContentsWord(); !reader.isContentsEnded(); val = reader.GetContentsWord() )
            {
              reader.ParseReal(val,Vector,Repeat,nnodes);
              for(int l = 0; l < Repeat; ++l)
              {
                for(int q = 0; q < (int)Vector.size(); ++q)
                {
                  xyz.push_back(Vector[q]);
                  if( xyz.size() == ndims )
                  {
                    new_nodes.push_back(CreateNode(xyz.data())->GetHandle());
                    xyz.clear();
                  }
                }
              }
            }
            reader.ReadCloseContents();
          }
          else 
          {
            reader.Report("Cannot find contents of XML tag");
            throw BadFile;
          }

          reader.CloseTag(TagNodes);
        }

        if( new_nodes.size() != nnodes )
        {
          reader.Report("Number of records for nodes %d do not match specified number of coordinates %d",new_nodes.size(),nnodes);
        }
      }

      { //Edges, Faces, Cells
        bool repeat = false;
          
        do
        {
          int nelems = 0, ntotconns = 0;
          bool matchelems = false;
          std::vector<HandleType> * elems;
          ElementType curtype;


          XMLReader::XMLTag TagElems;
          
          if( pass_tag )
          {
            TagElems = PassTag;
            pass_tag = false;
          }
          else
            TagElems = reader.OpenTag();

          if( !(TagElems.name == "Cells" || TagElems.name == "Faces" || TagElems.name == "Edges" ) )
          {
            reader.Report("Unexpected tag %s while waiting for either Cells or Faces or Edges",TagElems.name.c_str());
            throw BadFile;
          }
          if( TagElems.name != "Cells" ) repeat = true; else repeat = false;
          if( TagElems.name == "Cells" ) 
          {
            elems = &new_cells;
            curtype = CELL;
          }
          else if( TagElems.name == "Faces" ) 
          {
            elems = &new_faces;
            curtype = FACE;
          }
          else if( TagElems.name == "Edges" ) 
          {
            elems = &new_edges;
            curtype = EDGE;
          }

          for(int q = 0; q < TagElems.NumAttrib(); ++q)
          {
            XMLReader::XMLAttrib & attr = TagElems.GetAttib(q);
            if( attr.name == "Number" ) nelems = atoi(attr.value.c_str());
            else reader.Report("Unused attribute for %ss %s='%s'",TagElems.name.c_str(),attr.name.c_str(),attr.value.c_str());
          } 
          
          elems->reserve(nelems);

            
          HandleType * links[3] =
          {
            new_nodes.empty() ? NULL : &new_nodes[0],
            new_edges.empty() ? NULL : &new_edges[0],
            new_faces.empty() ? NULL : &new_faces[0]
          };

          
          XMLReader::XMLTag TagConns;
          for(TagConns = reader.OpenTag(); !TagConns.Finalize() && TagConns.name == "Connections"; reader.CloseTag(TagConns), TagConns = reader.OpenTag())
          {
            int nconns = 0;
            int offset = 0;
            ElementType subtype = NODE;
            for(int q = 0; q < TagConns.NumAttrib(); ++q)
            {
              XMLReader::XMLAttrib & attr = TagConns.GetAttib(q);
              if( attr.name == "Number" ) 
              {
                nconns = atoi(attr.value.c_str());
                matchelems = true;
              }
              else if( attr.name == "Type" ) subtype = reader.atoes(attr.value.c_str());
              else if( attr.name == "Offset" ) offset = atoi(attr.value.c_str());
              else reader.Report("Unused attribute for %ss %s='%s'",TagConns.name.c_str(),attr.name.c_str(),attr.value.c_str());
            }
            
            ntotconns += nconns;
            if( subtype >= curtype )
            {
              reader.Report("%ss cannot be constructed from %ss",ElementTypeName(curtype),ElementTypeName(subtype));
              throw BadFile;
            }
            reader.ReadOpenContents();
            ElementArray<Element> subarr(this);
            for(std::string val = reader.GetContentsWord(); !reader.isContentsEnded(); val = reader.GetContentsWord())
            {
              int num = atoi(val.c_str()), elem;
              subarr.clear();
              subarr.reserve(num);
              for(int q = 0; q < num && !reader.isContentsEnded(); ++q)
              {
                val = reader.GetContentsWord();
                elem = atoi(val.c_str());
                //subarr.push_back(ComposeHandle(subtype,elem-offset));
                subarr.push_back(links[ElementNum(subtype)][elem-offset]);
              }
              switch(curtype)
              {
              case EDGE:
                elems->push_back(CreateEdge(subarr.Convert<Node>()).first.GetHandle());
                break;
              case FACE:
                if( subtype == NODE )
                  elems->push_back(CreateFace(subarr.Convert<Node>()).first.GetHandle());
                else if( subtype == EDGE )
                  elems->push_back(CreateFace(subarr.Convert<Edge>()).first.GetHandle());
                break;
              case CELL:
                if( subtype == NODE )
                {
                  switch(subarr.size())
                  {
                  case 4: //Tetrahedron
                    {
                      const integer nodesnum[12] = {0,2,1,0,1,3,1,2,3,0,3,2};
										  const integer sizes[4] = {3,3,3,3};
                      elems->push_back(CreateCell(subarr.Convert<Node>(),nodesnum,sizes,4).first.GetHandle());
                    }
                    break;
                  case 5: //Pyramid
                    {
                      const integer nodesnum[16] = {0,4,3,0,1,4,1,2,4,3,4,2,0,3,2,1};
										  const integer sizes[5] = {3,3,3,3,4};
                      elems->push_back(CreateCell(subarr.Convert<Node>(),nodesnum,sizes,5).first.GetHandle());
                    }
                    break;
                  case 6: //Wedge or Prism
                    {
                      const integer nodesnum[18] = {0,2,5,3,1,4,5,2,0,3,4,1,3,5,4,0,1,2};
										  const integer sizes[5] = {4,4,4,3,3};
                      elems->push_back(CreateCell(subarr.Convert<Node>(),nodesnum,sizes,5).first.GetHandle());
                    }
                    break;
                  case 8: //Hexahedron
                    {
                      const integer nodesnum[24] = {0,4,7,3,1,2,6,5,0,1,5,4,3,7,6,2,0,3,2,1,4,5,6,7};
										  const integer sizes[6] = {4,4,4,4,4,4};
                      elems->push_back(CreateCell(subarr.Convert<Node>(),nodesnum,sizes,6).first.GetHandle());
                    }
                    break;
                  default: 
                    reader.Report("Sorry, no rule to convert %d nodes into a cell",subarr.size());
                    throw BadFile;
                    break;
                  }
                }
                else if( subtype == EDGE )
                {
                  reader.Report("Sorry, no rule to convert %d edges into a cell",subarr.size());
                  throw BadFile;
                }
                else if( subtype == FACE )
                  elems->push_back(CreateCell(subarr.Convert<Face>()).first.GetHandle());
                break;
              }
            }
            reader.ReadCloseContents();
          } 

          if( matchelems && nelems != ntotconns) reader.Report("Number %d of elements encountered do not match to the specified number %d",ntotconns,nelems);

          reader.CloseTag(TagElems);

          if( !TagConns.Finalize() )
          {
            PassTag = TagConns;
            pass_tag = true;
          }

        } while(repeat);
      }

      { //Sets and Data
        bool repeat = false;
        do
        {
          XMLReader::XMLTag TagSetsData;
          if( pass_tag )
          {
            TagSetsData = PassTag;
            pass_tag = false;
          }
          else TagSetsData = reader.OpenTag();
          
          if( TagSetsData.name == "Sets" )
          {
            repeat = true;
            HandleType * links[4] =
            {
              new_nodes.empty() ? NULL : &new_nodes[0],
              new_edges.empty() ? NULL : &new_edges[0],
              new_faces.empty() ? NULL : &new_faces[0],
              new_cells.empty() ? NULL : &new_cells[0]
            };
            int nsets = 0, sets_offset = 0;
            int nsets_read = 0;
            bool matchsets = false;
            for(int q = 0; q < TagSetsData.NumAttrib(); ++q)
            {
              XMLReader::XMLAttrib & attr = TagSetsData.GetAttib(q);
              if( attr.name == "Number" ) 
              {
                nsets = atoi(attr.value.c_str());
                matchsets = true;
              }
              else if( attr.name == "Offset" ) sets_offset = atoi(attr.value.c_str());
              else reader.Report("Unused attribute for %ss %s='%s'",TagSetsData.name.c_str(),attr.name.c_str(),attr.value.c_str());
            }
            new_sets.reserve(nsets);
            
            XMLReader::XMLTag Set;
            for(Set = reader.OpenTag(); !Set.Finalize() && Set.name == "Set"; reader.CloseTag(Set), Set = reader.OpenTag() )
            {
              int size = 0, offset = 0;
              std::string name;
              HandleType parent = InvalidHandle(), child = InvalidHandle(), sibling = InvalidHandle();
              ElementSet::ComparatorType comparator = ElementSet::UNSORTED_COMPARATOR;
              nsets_read++;
              for(int q = 0; q < Set.NumAttrib(); ++q)
              {
                XMLReader::XMLAttrib & attr = Set.GetAttib(q);
                if( attr.name == "Size" ) size = atoi(attr.value.c_str());
                else if( attr.name == "Offset" ) offset = atoi(attr.value.c_str());
                else if( attr.name == "Name" ) name = attr.value;
                else if( attr.name == "Parent" ) 
                {
                  if( attr.value != "Unset" )
                    parent = ComposeHandle(ESET,atoi(attr.value.c_str()));
                }
                else if( attr.name == "Sibling" ) 
                {
                  if( attr.value != "Unset" )
                    sibling = ComposeHandle(ESET,atoi(attr.value.c_str()));
                }
                else if( attr.name == "Child" ) 
                {
                  if( attr.value != "Unset" )
                    child = ComposeHandle(ESET,atoi(attr.value.c_str()));
                }
                else if( attr.name == "Comparator" )
                {
                  if( attr.value == "Unsorted" ) comparator = ElementSet::UNSORTED_COMPARATOR;
                  else if( attr.value == "Identificator" ) comparator = ElementSet::GLOBALID_COMPARATOR;
                  else if( attr.value == "Centroid" ) comparator = ElementSet::CENTROID_COMPARATOR;
                  else if( attr.value == "Hierarchy" ) comparator = ElementSet::HIERARCHY_COMPARATOR;
                  else if( attr.value == "Handle" ) comparator = ElementSet::HANDLE_COMPARATOR;
                  else reader.Report("Unexpected comparator type %s for attribute Comparator, expected Unsorted,Identificator,Centroid,Hierarchy,Handle",attr.value.c_str());
                }
                else reader.Report("Unused attribute for %ss %s='%s'",Set.name.c_str(),attr.name.c_str(),attr.value.c_str());
              } 
              ElementSet s = CreateSet(name).first;
              new_sets.push_back(s.GetHandle());
              set_comparators.push_back(comparator);
              Element::adj_type & hc = HighConn(s.GetHandle());
              Element::adj_type & lc = LowConn(s.GetHandle());
              //here we believe that all the connections are consistent
              //hc.resize(ElementSet::high_conn_reserved);
              if( parent != InvalidHandle() ) ElementSet::hParent(hc) = parent;
              if( child != InvalidHandle() ) ElementSet::hChild(hc) = child;
              if( sibling != InvalidHandle() ) ElementSet::hSibling(hc) = sibling;
              //s.SortSet(comparator);
              reader.ReadOpenContents();
              for(std::string val = reader.GetContentsWord(); !reader.isContentsEnded(); val = reader.GetContentsWord())
              {
                std::vector<std::pair<ElementType,int> > Vector;
                int Repeat;
                reader.ParseReference(val,Vector,Repeat,size);
                if( !Vector.empty() ) for(int l = 0; l < Repeat; ++l) 
                {
                  for(int q = 0; q < (int)Vector.size(); ++q)
                  {
                    if( ElementNum(Vector[q].first) < 4 )
                      lc.push_back(links[ElementNum(Vector[q].first)][Vector[q].second-offset]);
                    else if( ElementNum(Vector[q].first) == 4 ) //Set
                      lc.push_back(ComposeHandle(Vector[q].first,Vector[q].second-offset));
                    else //Mesh
                      lc.push_back(GetHandle());
                  }
                    //s.PutElement(links[ElementNum(Vector[q].first)][Vector[q].second-offset]);
                  //s.PutElements(Vector[0],(enumerator)Vector.size());
                }
              }
              reader.ReadCloseContents();
            }
            
            reader.CloseTag(TagSetsData);
            if( !Set.Finalize() )
            {
              PassTag = Set;
              pass_tag = true;
            }
            
            if( matchsets && nsets != nsets_read ) reader.Report("Number %d of XML tags Set read do not match to the number %d specified",nsets_read,nsets);

            //correct links between sets
            for(int q = 0; q < (int)new_sets.size(); ++q)
            {
              Element::adj_type & lc = HighConn(new_sets[q]);
              for(Element::adj_type::iterator jt = lc.begin(); jt != lc.end(); ++jt)
              {
                if( GetHandleElementType(*jt) == ESET )
                  *jt = new_sets[GetHandleID(*jt)];
              }
              Element::adj_type & hc = HighConn(new_sets[q]);
              if( ElementSet::hParent(hc) != InvalidHandle() ) ElementSet::hParent(hc) = new_sets[GetHandleID(ElementSet::hParent(hc))-sets_offset];
              if( ElementSet::hChild(hc) != InvalidHandle() ) ElementSet::hChild(hc) = new_sets[GetHandleID(ElementSet::hChild(hc))-sets_offset];
              if( ElementSet::hSibling(hc) != InvalidHandle() ) ElementSet::hSibling(hc) = new_sets[GetHandleID(ElementSet::hSibling(hc))-sets_offset];
              if( set_comparators[q] != ElementSet::UNSORTED_COMPARATOR ) ElementSet(this,new_sets[q]).SortSet(set_comparators[q]);
            }
          }
          else if( TagSetsData.name == "Data" )
          {
            repeat = false;
            HandleType * links[5] =
            {
              new_nodes.empty() ? NULL : &new_nodes[0],
              new_edges.empty() ? NULL : &new_edges[0],
              new_faces.empty() ? NULL : &new_faces[0],
              new_cells.empty() ? NULL : &new_cells[0],
              new_sets.empty() ? NULL : &new_sets[0]
            };
            int ndata = 0, ndatapass = 0;
            bool matchndata = false;

            for(int q = 0; q < TagSetsData.NumAttrib(); ++q)
            {
              XMLReader::XMLAttrib & attr = TagSetsData.GetAttib(q);
              if( attr.name == "Number" ) 
              {
                ndata = atoi(attr.value.c_str());
                matchndata = true;
              }
              else reader.Report("Unused attribute for %ss %s='%s'",TagSetsData.name.c_str(),attr.name.c_str(),attr.value.c_str());
            }

            XMLReader::XMLTag TagDataSet;
            for(TagDataSet = reader.OpenTag(); !TagDataSet.Finalize() && TagDataSet.name == "DataSet"; reader.CloseTag(TagDataSet), TagDataSet = reader.OpenTag())
            {
              ++ndatapass;
              int sparse_read = 2, offset = 0;
              std::string tagname = "", setname = "", meshname = "";
              Mesh * remote_mesh = NULL;
              ElementType etype = NONE;
              for(int q = 0; q < TagDataSet.NumAttrib(); ++q)
              {
                XMLReader::XMLAttrib & attr = TagDataSet.GetAttib(q);
                if( attr.name == "SetType" ) 
                {
                  if( attr.value == "Cells" ) etype = CELL;
                  else if( attr.value == "Faces" ) etype = FACE;
                  else if( attr.value == "Edges" ) etype = EDGE;
                  else if( attr.value == "Nodes" ) etype = NODE;
                  else if( attr.value == "Sets" ) etype = ESET;
                  else if( attr.value == "Mesh" ) etype = MESH;
                  else if( attr.value == "SetData" ) etype = NONE;
                }
                else if( attr.name == "TagName" ) tagname = attr.value;
                else if( attr.name == "SetName" ) setname = attr.value;
                else if( attr.name == "MeshName" ) meshname = attr.value;
                else if( attr.name == "Sparse" ) sparse_read = reader.ParseBool(attr.value);
                else if( attr.name == "Offset" ) offset = atoi(attr.value.c_str());
                else reader.Report("Unused attribute for %ss %s='%s'",TagDataSet.name.c_str(),attr.name.c_str(),attr.value.c_str());
              }
              
              if( tagname == "" )
              {
                reader.Report("DataSet had no attribute TagName");
                throw BadFile;
              }
              if( !HaveTag(tagname) )
              {
                reader.Report("Tag %s do not exist",tagname.c_str());
                throw BadFile;
              }
              if( meshname != "" )
              {
                remote_mesh = Mesh::GetMesh(meshname);
                if( !remote_mesh )
                  reader.Report("Remote mesh %s do not exist, you should create it first inside of your application",meshname.c_str());
              }

              Tag t = GetTag(tagname);
              HandleType * set_elems = NULL, *it = NULL;
              enumerator set_size = 0;

              if( setname != "" )
              {
                if( etype != NONE ) reader.Report("Warning: SetType should be SetData for data specified for set.");
                ElementSet s = GetSet(setname);
                if( !s.isValid() )
                {
                  reader.Report("Set %s is not defined on the mesh",setname.c_str());
                  throw BadFile;
                }
                set_elems = s.getHandles();
                set_size = s.Size(); 
                if( sparse_read == 2 ) sparse_read = 0;
              }
              else
              {
                if( etype == NONE )
                {
                  for(ElementType test = NODE; test <= MESH; test = NextElementType(test) )
                    if( t.isDefined(test) ) etype |= test;
                  if( !OneType(etype) )
                  {
                    reader.Report("You must explicitly specify type of elements for which"
                                  "the data is written since tag %s is defined on multiple types of elements",tagname.c_str());
                    throw BadFile;
                  }
                }
                if( sparse_read == 2 )
                {
                  if( t.isSparse(etype) ) sparse_read = 1; //Expect data to be listed in sparse manner for sparse tag
                  else sparse_read = 0;
                }
                switch(etype)
                {
                case NODE: set_elems = &new_nodes[0]; set_size = (enumerator)new_nodes.size(); break;
                case EDGE: set_elems = &new_edges[0]; set_size = (enumerator)new_edges.size(); break;
                case FACE: set_elems = &new_faces[0]; set_size = (enumerator)new_faces.size(); break;
                case CELL: set_elems = &new_cells[0]; set_size = (enumerator)new_cells.size(); break;
                case ESET: set_elems = &new_sets[0]; set_size = (enumerator)new_sets.size(); break;
                case MESH: set_elems = &handle; set_size = 1; break;
                }
              }
              it = set_elems;
              reader.ReadOpenContents();
              for(std::string val = reader.GetContentsWord(); !reader.isContentsEnded(); val = reader.GetContentsWord())
              {
                if( sparse_read )
                {
                  it = set_elems + atoi(val.c_str());
                  val = reader.GetContentsWord();
                }
                switch(t.GetDataType())
                {
                case DATA_REAL:
                  {
                    Storage::real_array data = RealArray(*it,t);
                    std::vector<Storage::real> Vector; int Repeat;
                    reader.ParseReal(val,Vector,Repeat,set_size);
                    if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                    else if( t.GetSize() != Vector.size()*Repeat && sparse_read )
                    {
                      reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                      throw BadFile;
                    }
                      
                    for(int l = 0; l < Repeat; ++l)
                    {
                      for(int q = 0; q < (int)Vector.size(); ++q)
                      {
                        data[(q + l*((int)Vector.size()))%data.size()] = Vector[q];
                        if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                        {
                          ++it;
                          if( ((int)(it-set_elems)) < (int)set_size ) data = RealArray(*it,t);
                        }
                      }
                    }
                  }
                  break;
                case DATA_INTEGER:
                  {
                    Storage::integer_array data = IntegerArray(*it,t);
                    std::vector<Storage::integer> Vector; int Repeat;
                    reader.ParseInteger(val,Vector,Repeat,set_size);
                    if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                    else if( t.GetSize() != Vector.size()*Repeat && sparse_read)
                    {
                      reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                      throw BadFile;
                    }
                    for(int l = 0; l < Repeat; ++l)
                    {
                      for(int q = 0; q < (int)Vector.size(); ++q)
                      {
                        data[(q + l*((int)Vector.size()))%data.size()] = Vector[q];
                        if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                        {
                          ++it;
                          if( ((int)(it-set_elems)) < (int)set_size ) data = IntegerArray(*it,t);
                        }
                      }
                    }
                  }
                  break;
                case DATA_BULK:
                  {
                    Storage::bulk_array data = BulkArray(*it,t);
                    std::string Vector; int Repeat;
                    reader.ParseBulk(val,Vector,Repeat,set_size);
                    if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                    else if( t.GetSize() != Vector.size()*Repeat && sparse_read)
                    {
                      reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                      throw BadFile;
                    }
                    for(int l = 0; l < Repeat; ++l)
                    {
                      for(int q = 0; q < (int)Vector.size(); ++q)
                      {
                        data[(q + l*((int)Vector.size()))%data.size()] = Vector[q];
                        if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                        {
                          ++it;
                          if( ((int)(it-set_elems)) < (int)set_size ) data = BulkArray(*it,t);
                        }
                      }
                    }
                  }
                  break;
                case DATA_REFERENCE:
                  {
                    Storage::reference_array data = ReferenceArray(*it,t);
                    std::vector<std::pair<ElementType,int> > Vector; int Repeat;
                    reader.ParseReference(val,Vector,Repeat,set_size);
                    if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                    else if( t.GetSize() != Vector.size()*Repeat && sparse_read)
                    {
                      reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                      throw BadFile;
                    }
                    for(int l = 0; l < Repeat; ++l)
                    {
                      for(int q = 0; q < (int)Vector.size(); ++q)
                      {
                        data[(q + l*((int)Vector.size()))%data.size()] = Element(this,links[ElementNum(Vector[q].first)][Vector[q].second-offset]);
                        if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                        {
                          ++it;
                          if( ((int)(it-set_elems)) < (int)set_size ) data = ReferenceArray(*it,t);
                        }
                      }
                    }
                  }
                  break;
                case DATA_REMOTE_REFERENCE:
                  {
                    Storage::remote_reference_array data = RemoteReferenceArray(*it,t);
                    if( remote_mesh != NULL )
                    {
                      std::vector<std::pair<ElementType,int> > Vector; int Repeat;
                      reader.ParseReference(val,Vector,Repeat,set_size);
                      if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                      else if( t.GetSize() != Vector.size()*Repeat && sparse_read )
                      {
                        reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                        throw BadFile;
                      }
                      for(int l = 0; l < Repeat; ++l)
                      {
                        for(int q = 0; q < (int)Vector.size(); ++q)
                        {
                          data.at((q + l*((int)Vector.size()))%data.size()).first  = remote_mesh;
                          data.at((q + l*((int)Vector.size()))%data.size()).second = ComposeHandle(Vector[q].first,Vector[q].second-offset);
                          if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                          {
                            ++it;
                            if( ((int)(it-set_elems)) < (int)set_size ) data = RemoteReferenceArray(*it,t);
                          }
                        }
                      }
                    }
                    else
                    {
                      std::vector<std::pair<std::string,std::pair<ElementType,int> > > Vector; int Repeat;
                      reader.ParseRemoteReference(val,Vector,Repeat,set_size);
                      if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                      else if( t.GetSize() != Vector.size()*Repeat && sparse_read )
                      {
                        reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                        throw BadFile;
                      }
                      for(int l = 0; l < Repeat; ++l)
                      {
                        for(int q = 0; q < (int)Vector.size(); ++q)
                        {
                          Mesh * m = GetMesh(Vector[q].first);
                          if( m == NULL ) reader.Report("Cannot find remote mesh %s, you should create it first inside of your application",Vector[q].first.c_str());
                          data.at((q + l*((int)Vector.size()))%data.size()).first = m;
                          data.at((q + l*((int)Vector.size()))%data.size()).second = ComposeHandle(Vector[q].second.first,Vector[q].second.second-offset);
                          if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                          {
                            ++it;
                            if( ((int)(it-set_elems)) < (int)set_size ) data = RemoteReferenceArray(*it,t);
                          }
                        }
                      }
                    }
                  }
                  break;
#if defined(USE_AUTODIFF)
                case DATA_VARIABLE:
                  {
                    Storage::var_array data = VariableArray(*it,t);
                    std::vector<Storage::var> Vector; int Repeat;
                    reader.ParseVariable(val,Vector,Repeat,set_size);
                    if( t.GetSize() == ENUMUNDEF ) data.resize((enumerator)Vector.size()*Repeat);
                    else if( t.GetSize() != Vector.size()*Repeat && sparse_read )
                    {
                      reader.Report("Cannot write record of size %d into tag data of size %d",Vector.size()*Repeat,t.GetSize());
                      throw BadFile;
                    }
                    for(int l = 0; l < Repeat; ++l)
                    {
                      for(int q = 0; q < (int)Vector.size(); ++q)
                      {
                        data[(q + l*((int)Vector.size()))%data.size()] = Vector[q];
                        if( !sparse_read && t.GetSize() != ENUMUNDEF && (q + l*((int)Vector.size())+1)%t.GetSize() == 0 )
                        {
                          ++it;
                          if( ((int)(it-set_elems)) < (int)set_size ) data = VariableArray(*it,t);
                        }
                      }
                    }
                  }
                  break;
#endif
                }
              } 
              reader.ReadCloseContents();
            }

            if( matchndata && ndata != ndatapass ) reader.Report("warning: Number %d of DataSet tags encountered do not match to the specified %d",ndatapass,ndata);

            reader.CloseTag(TagSetsData);

            if( !TagDataSet.Finalize() )
            {
              PassTag = TagDataSet;
              pass_tag = true;
            }
          }
          else
          {
            reader.Report("Unexpected tag %s, expected Sets or Data",TagSetsData.name.c_str());
            throw BadFile;
          }
        } while(repeat);
      }
      RepairGeometricTags();

      if( repair_orientation )
      {
        int numfixed = 0;
        for(int q = 0; q < (int)new_faces.size(); ++q)
          if( Face(this,new_faces[q]).FixNormalOrientation() ) numfixed++;
        if( numfixed ) std::cout << "Fixed orientation of " << numfixed << " faces" << std::endl;
      }
    }
    reader.CloseTag(TagParallelMesh);
    infile.close();
  }

  void Mesh::SaveXML(std::string File)
  {
    std::fstream fout(File.c_str(),std::ios::out);
    fout << "<ParallelMesh>\n";
    fout << "\t<Mesh>\n";
    fout << "\t\t<Tags>\n";
    for(int k = 0; k < (int)tags.size(); ++k)
    {
      if( tags[k].GetTagName().substr(0,9) == "PROTECTED" ) continue;
      std::string names[6] = {"Nodes","Edges","Faces","Cells","Sets","Mesh"};
      std::string definition = "", sparse = "", type = "";
      switch(tags[k].GetDataType())
      {
      case DATA_REAL: type = "Real"; break;
      case DATA_INTEGER: type = "Integer"; break;
      case DATA_BULK: type = "Bulk"; break;
      case DATA_REFERENCE: type = "Reference"; break;
      case DATA_REMOTE_REFERENCE: type = "RemoteReference"; break;
#if defined(USE_AUTODIFF)
      case DATA_VARIABLE: type = "Variable"; break;
#endif
      }
      for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype) )
      {
        if( tags[k].isDefined(etype) ) definition += names[ElementNum(etype)] + ",";
        if( tags[k].isSparse(etype) ) sparse += names[ElementNum(etype)] + ",";
      }
      if( !definition.empty() ) definition.resize(definition.size()-1); //remove trailing comma
      if( !sparse.empty() ) sparse.resize(sparse.size()-1); //remove trailing comma
      if( sparse == "" ) sparse = "None";
      fout << "\t\t\t<Tag Name  =\"" << tags[k].GetTagName() << "\"\n";
      if( tags[k].GetSize() !=ENUMUNDEF )
      fout << "\t\t\t     Size  =\"" << tags[k].GetSize() << "\"\n";
      fout << "\t\t\t     Type  =\"" << type << "\"\n";
      fout << "\t\t\t     Sparse=\"" << sparse << "\"\n";
      fout << "\t\t\t     Definition=\"" << definition << "\"/>\n";
    }
    Tag idx = CreateTag("TEMPORARY_XML_ENUMERATOR",DATA_INTEGER,MESH|ESET|CELL|FACE|EDGE|NODE,NONE,1);
    Integer(GetHandle(),idx) = 0;
    fout << "\t\t</Tags>\n";
    fout << "\t\t<Nodes Number=\"" << NumberOfNodes() << "\" Dimensions=\"" << GetDimensions() << "\">\n";
    fout << "\t\t\t<![CDATA[\n";
    dynarray<Storage::real,3> xyz(GetDimensions());
    int cnt = 0;
    for(iteratorNode it = BeginNode(); it != EndNode(); ++it)
    {
      fout << "\t\t\t";
      it->Centroid(xyz.data());
      for(int q = 0; q < GetDimensions(); ++q)
        fout << xyz[q] << " ";
      fout << "\n";
      it->Integer(idx) = cnt++;
    }
    fout << "\t\t\t]]>\n";
    fout << "\t\t</Nodes>\n";
    fout << "\t\t<Edges>\n";
    fout << "\t\t\t<Connections Number=\"" << NumberOfEdges() << "\" Type=\"Nodes\">\n";
    fout << "\t\t\t<![CDATA[\n";
    cnt = 0;
    for(iteratorEdge it = BeginEdge(); it != EndEdge(); ++it)
    {
      ElementArray<Node> nodes = it->getNodes();
      fout << "\t\t\t" << nodes.size();
      for(ElementArray<Node>::iterator jt = nodes.begin(); jt != nodes.end(); ++jt)
        fout << " " << jt->Integer(idx);
      fout << "\n";
      it->Integer(idx) = cnt++;
    }
    fout << "\t\t\t]]>\n";
    fout << "\t\t\t</Connections>\n";
    fout << "\t\t</Edges>\n";
    fout << "\t\t<Faces>\n";
    fout << "\t\t\t<Connections Number=\"" << NumberOfFaces() << "\" Type=\"Edges\">\n";
    fout << "\t\t\t<![CDATA[\n";
    cnt = 0;
    for(iteratorFace it = BeginFace(); it != EndFace(); ++it)
    {
      ElementArray<Edge> edges = it->getEdges();
      fout << "\t\t\t" << edges.size();
      for(ElementArray<Edge>::iterator jt = edges.begin(); jt != edges.end(); ++jt)
        fout << " " << jt->Integer(idx);
      fout << "\n";
      it->Integer(idx) = cnt++;
    }
    fout << "\t\t\t]]>\n";
    fout << "\t\t\t</Connections>\n";
    fout << "\t\t</Faces>\n";
    fout << "\t\t<Cells>\n";
    fout << "\t\t\t<Connections Number=\"" << NumberOfFaces() << "\" Type=\"Faces\">\n";
    fout << "\t\t\t<![CDATA[\n";
    cnt = 0;
    for(iteratorCell it = BeginCell(); it != EndCell(); ++it)
    {
      ElementArray<Face> faces = it->getFaces();
      fout << "\t\t\t" << faces.size();
      for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
        fout << " " << jt->Integer(idx);
      fout << "\n";
      it->Integer(idx) = cnt++;
    }
    fout << "\t\t\t]]>\n";
    fout << "\t\t\t</Connections>\n";
    fout << "\t\t</Cells>\n";
    fout << "\t\t<Sets>\n";
    cnt = 0;
    for(iteratorSet it = BeginSet(); it != EndSet(); ++it)
      it->Integer(idx) = cnt++;
    for(iteratorSet it = BeginSet(); it != EndSet(); ++it)
    {
      fout << "\t\t\t<Set Name=\"" << it->GetName() << "\"\n";
      if( it->HaveParent() )
      fout << "\t\t\t     Parent=\"" << it->GetParent()->Integer(idx) << "\"\n";
      if( it->HaveChild() )
      fout << "\t\t\t     Child=\"" << it->GetChild()->Integer(idx) << "\"\n";
      if( it->HaveSibling() )
      fout << "\t\t\t     Sibling=\"" << it->GetSibling()->Integer(idx) << "\"\n";
      if( it->GetComparator() != ElementSet::UNSORTED_COMPARATOR )
      {
        fout << "\t\t\t     Comparator=\"";
        switch(it->GetComparator())
        {
        case ElementSet::GLOBALID_COMPARATOR: fout << "Identificator"; break;
        case ElementSet::CENTROID_COMPARATOR: fout << "Centroid"; break;
        case ElementSet::HIERARCHY_COMPARATOR: fout << "Hierarchy"; break;
        case ElementSet::HANDLE_COMPARATOR: fout << "Handle"; break;
        }
        fout << "\"\n";
      }
      fout << "\t\t\t     Size=\"" << it->Size() << "\">\n";
      fout << "\t\t\t<![CDATA[\n";
      int endl_count = 0;
      fout << "\t\t\t";
      for(ElementSet::iterator jt = it->Begin(); jt != it->End(); ++jt)
      {
        switch(jt->GetElementType())
        {
        case NODE: fout << "Node:"; break;
        case EDGE: fout << "Edge:"; break;
        case FACE: fout << "Face:"; break;
        case CELL: fout << "Cell:"; break;
        case ESET: fout << "Set:"; break;
        case MESH: fout << "Mesh:"; break;
        }
        fout << jt->Integer(idx) << " ";
        endl_count++;
        if( endl_count % 8 == 0 )
          fout << "\n\t\t\t";
      }
      fout << "\n\t\t\t]]>\n";
      fout << "\t\t\t</Set>\n";
    }
    fout << "\t\t</Sets>\n";
    fout << "\t\t<Data>\n";
    for(iteratorTag t = BeginTag(); t != EndTag(); ++t)
    {
      if( *t == idx ) continue;
      if( t->GetTagName().substr(0,9) == "PROTECTED" ) continue;
      for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype)) if( t->isDefined(etype) )
      {
        fout << "\t\t\t<DataSet TagName=\"" << t->GetTagName() << "\"\n";
        fout << "\t\t\t         SetType=\"";
        switch(etype)
        {
        case NODE: fout << "Nodes"; break;
        case EDGE: fout << "Edges"; break;
        case FACE: fout << "Faces"; break;
        case CELL: fout << "Cells"; break;
        case ESET: fout << "Sets"; break;
        case MESH: fout << "Mesh"; break;
        }
        fout << "\">\n";
        fout << "\t\t\t<![CDATA[\n\t\t\t";
        int endl_count = 0;
        switch(t->GetDataType())
        {
        case DATA_REAL:
          for(iteratorStorage jt = Begin(etype); jt != End(); ++jt) if( jt->HaveData(*t) )
          {
            Storage::real_array data = jt->RealArray(*t);
            if( !data.empty() )
            {
              if( t->isSparse(etype) ) 
              {
                fout << jt->Integer(idx) << " ";
                endl_count++;
              }
              if( data.size() == 1 ) 
              {
                fout << data[0];
              }
              else
              {
                fout << "{";
                for(int q = 0; q < (int)data.size()-1; ++q)
                  fout << data[q] << ",";
                fout << data.back();
                fout << "}";
              }
              endl_count+=data.size();
              if( endl_count > 8 )
              {
                fout << "\n\t\t\t";
                endl_count = 0;
              }
              else fout << " ";
            }
          }
          break;
        case DATA_INTEGER:
          for(iteratorStorage jt = Begin(etype); jt != End(); ++jt) if( jt->HaveData(*t) )
          {
            Storage::integer_array data = jt->IntegerArray(*t);
            if( !data.empty() )
            {
              if( t->isSparse(etype) ) 
              {
                fout << jt->Integer(idx) << " ";
                endl_count++;
              }
              if( data.size() == 1 ) fout << data[0];
              else
              {
                fout << "{";
                for(int q = 0; q < (int)data.size()-1; ++q)
                  fout << data[q] << ",";
                fout << data.back();
                fout << "}";
              }
              endl_count+=data.size();
              if( endl_count > 8 )
              {
                fout << "\n\t\t\t";
                endl_count = 0;
              }
              else fout << " ";
              //fout << std::endl;
            }
          }
          break;
        case DATA_BULK:
          for(iteratorStorage jt = Begin(etype); jt != End(); ++jt) if( jt->HaveData(*t) )
          {
            Storage::bulk_array data = jt->BulkArray(*t);
            if( !data.empty() )
            {
              if( t->isSparse(etype) ) 
              {
                fout << jt->Integer(idx) << " ";
                endl_count++;
              }
              if( data.size() == 1 ) fout << CharToHex(data[0]);
              else
              {
                fout << "{";
                for(int q = 0; q < (int)data.size()-1; ++q)
                  fout << CharToHex(data[q]) << ",";
                fout << CharToHex(data.back());
                fout << "}";
              }
              endl_count+=data.size();
              if( endl_count > 8 )
              {
                fout << "\n\t\t\t";
                endl_count = 0;
              }
              else fout << " ";
              //fout << std::endl;
            }
          }
          break;
        case DATA_REFERENCE:
          for(iteratorStorage jt = Begin(etype); jt != End(); ++jt) if( jt->HaveData(*t) )
          {
            Storage::reference_array data = jt->ReferenceArray(*t);
            if( !data.empty() )
            {
              if( t->isSparse(etype) ) 
              {
                fout << jt->Integer(idx) << " ";
                endl_count++;
              }
              if( data.size() == 1 ) 
              {
                if( data[0].isValid() )
                  fout << ReferenceToString(data.at(0),data[0].Integer(idx));
                else
                  fout << "None:0";
              }
              else
              {
                fout << "{";
                for(int q = 0; q < (int)data.size()-1; ++q)
                {
                  if( data[q].isValid() )
                    fout << ReferenceToString(data.at(q),data[q].Integer(idx)) << ",";
                  else
                    fout << "None:0,";
                }
                if( data[data.size()-1].isValid() )
                  fout << ReferenceToString(data.at(data.size()-1),data[data.size()-1].Integer(idx));
                else
                  fout << "None:0";
                fout << "}";
              }
              endl_count+=data.size();
              if( endl_count > 8 )
              {
                fout << "\n\t\t\t";
                endl_count = 0;
              }
              else fout << " ";
              //fout << std::endl;
            }
          }
          break;
        case DATA_REMOTE_REFERENCE:
          for(iteratorStorage jt = Begin(etype); jt != End(); ++jt) if( jt->HaveData(*t) )
          {
            Storage::remote_reference_array data = jt->RemoteReferenceArray(*t);
            if( !data.empty() )
            {
              if( t->isSparse(etype) ) 
              {
                fout << jt->Integer(idx) << " ";
                endl_count++;
              }
              if( data.size() == 1 ) 
              {
                if( data[0].isValid() )
                  fout << data.at(0).first->GetMeshName() << ":" << ReferenceToString(data.at(0).second,data[0].Integer(idx));
                else
                  fout << ":None:0";
              }
              else
              {
                fout << "{";
                for(int q = 0; q < (int)data.size()-1; ++q)
                {
                  if( data[q].isValid() )
                    fout << data.at(q).first->GetMeshName() << ":" << ReferenceToString(data.at(q).second,data[q].Integer(idx)) << ",";
                  else
                    fout << ":None:0,";
                }
                if( data[data.size()-1].isValid() )
                  fout << data.at(data.size()-1).first->GetMeshName() << ":" << ReferenceToString(data.at(data.size()-1).second,data[data.size()-1].Integer(idx));
                else
                  fout << ":None:0";
                fout << "}";
              }
              endl_count+=data.size();
              if( endl_count > 8 )
              {
                fout << "\n\t\t\t";
                endl_count = 0;
              }
              else fout << " ";
              //fout << std::endl;
            }
          }
          break;
#if defined(USE_AUTODIFF)
        case DATA_VARIABLE:
          for(iteratorStorage jt = Begin(etype); jt != End(); ++jt) if( jt->HaveData(*t) )
          {
            Storage::var_array data = jt->VariableArray(*t);
            if( !data.empty() )
            {
              if( t->isSparse(etype) ) 
              {
                fout << jt->Integer(idx) << " ";
                endl_count++;
              }
              if( data.size() == 1 ) fout << VariableToString(data.at(0));
              else
              {
                fout << "{";
                for(int q = 0; q < (int)data.size()-1; ++q)
                  fout << VariableToString(data.at(q)) << ",";
                fout << VariableToString(data.back());
                fout << "}";
              }
              endl_count+=data.size();
              if( endl_count > 8 )
              {
                fout << "\n\t\t\t";
                endl_count = 0;
              }
              else fout << " ";
              //fout << std::endl;
            }
          }
          break;
#endif
        }
        fout << "\n\t\t\t]]>\n";
        fout << "\t\t\t</DataSet>\n";
      }
    }
    DeleteTag(idx);
    fout << "\t\t</Data>\n";
    fout << "\t</Mesh>\n";
    fout << "</ParallelMesh>\n";
    fout.close();
  }
}

#endif