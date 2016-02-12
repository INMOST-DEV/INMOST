#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "inmost.h"
#include <stdarg.h> //for va_list

namespace INMOST
{

  static int get_priority(char c)
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

  static int ConvertHex(char in)
  {
    int ret = tolower(in) - 48;
    if( ret > 10 ) ret -= 7;
    if( ret > 15 ) ret -= 32;
    return ret;
  }
  static char atoc(const char * str)
  {
    return (char)(ConvertHex(str[0])*16 + ConvertHex(str[1]));
  }

  std::string CharToHex(char c)
  {
    char const hex[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B','C','D','E','F'};
    std::string str = "";
    str.append(&hex[(c & 0xF0) >> 4], 1);
    str.append(&hex[c & 0xF], 1);
    return str;
  }

#if defined(USE_MESH)
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
#endif

  static char * sstrip(char * str)
  {
	  int text_start = 0, text_end = (int)strlen(str);
	  for(text_start = 0; isspace(str[text_start]) && text_start < text_end; text_start++);
	  if( text_start == text_end ) return str+text_start;
	  for(text_end = text_end-1; isspace(str[text_end]) && text_end > text_start; text_end--);
	  str[text_end+1] = '\0';
	  return str+text_start;
  }

  static std::string sstrip(const std::string & input)
  {

	  char temp[2048];
	  strcpy(temp,input.c_str());
	  return std::string(sstrip(temp));
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

  std::vector<std::string> XMLReader::Interpreter::Expand(const std::string & input) const
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

  

  std::vector<std::string> XMLReader::Interpreter::MakePolish(const std::vector<std::string> & input) 
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
  void XMLReader::Interpreter::Print(const std::vector<std::string> & polish) const
  {
    for(int k = 0; k < (int)polish.size(); ++k)
      std::cout << polish[k] << " ";
    std::cout << std::endl;
  }
  double XMLReader::Interpreter::Run(const std::vector<std::string> & polish)
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

  XMLReader::Interpreter::Interpreter() :error_state(false) {}
  XMLReader::Interpreter::Interpreter(const Interpreter & b) : error_state(b.error_state) {}
  XMLReader::Interpreter & XMLReader::Interpreter::operator = (Interpreter const & b)
  {
    error_state = b.error_state;
    return * this;
  }

  double XMLReader::Interpreter::Evaluate(const std::string & str)
  {
    //const char * debug_str = str.c_str();
    std::vector<std::string> decompose = Expand(str);
    std::vector<std::string> polish = MakePolish(decompose);
    //Print(polish);
    return Run(polish);
  }
  bool XMLReader::Interpreter::isError() {return error_state;}
  void XMLReader::Interpreter::ClearError() {error_state = false;}

  XMLReader::Stream & XMLReader::get_Stream() {return inp.back();}
  const XMLReader::Stream & XMLReader::get_Stream() const {return inp.back();}
  std::istream & XMLReader::get_iStream() {return *inp.back().s;}
  const std::istream & XMLReader::get_iStream() const {return *inp.back().s;}

  XMLReader::XMLReader(const XMLReader & other) {}
  XMLReader & XMLReader::operator =(XMLReader & other) {return *this;}

  char XMLReader::GetChar()
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

  void XMLReader::RetChar()
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
  void XMLReader::SkipComments(State RetState)
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

  std::string XMLReader::StateName(State s) const
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


  void XMLReader::Report(const char * fmt, ...) const
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

  XMLReader::XMLReader(std::string sourcename, std::istream & input) :intrp(),_state(Intro)
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

  void XMLReader::PushStream(std::string file) 
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
  void XMLReader::PopStream() 
  {
    if( inp.size() > 1 )
      delete static_cast<std::fstream *>(inp.back().s);
    inp.pop_back();
    if( _state == EndOfFile && !inp.empty() )
      _state = Intro;
  }

  std::string XMLReader::ReadOpenTag()
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
  int XMLReader::ReadCloseTag()
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
  bool XMLReader::isTagFinish() const {return _state == ReadCloseTagSlash;}

  bool XMLReader::ReadFinishTag(std::string TagName)
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

  std::string XMLReader::AttributeName()
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

  std::string XMLReader::AttributeValue()
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

  bool XMLReader::isTagEnded() const {return _state == EndTag;}

  bool XMLReader::ReadOpenContents()
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

  std::string XMLReader::GetContentsWord()
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
  bool XMLReader::ReadCloseContents()
  {
    if( _state == EndContents && GetChar() == '>' )
    {
      _state = Intro;
      return true;
    }
    return false;
  }

  bool XMLReader::isContentsEnded() const
  {
    return _state == EndContents;
  }

  bool XMLReader::isFailure() const {return _state == Failure;}

  bool XMLReader::isEof() const {return _state == EndOfFile;}

#if defined(USE_MESH)
  INMOST::ElementType XMLReader::atoes(const char * _str)
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

  INMOST::ElementType XMLReader::atoe(const char * _str)
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

  std::pair<INMOST::ElementType,int> XMLReader::atoh(const char * _str)
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

  std::pair<std::string,std::pair<INMOST::ElementType,int> > XMLReader::atorh(const char * _str)
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
#endif

#if defined(USE_AUTODIFF)
  INMOST::Storage::var XMLReader::atov(const char * _str)
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

  int XMLReader::EvaluateExpression(std::string expression)
  {
    //port short RPN version later from DiscretizationToolkit/GPRS/interpreter.cpp
    return (int)intrp.Evaluate(expression);
    return atoi(expression.c_str()); //stupid version for now
  }

  int XMLReader::ConvertMultiplier(std::string expression, int SetSize)
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

  void XMLReader::SplitValueMultiplier(std::string expression, std::string & value, std::string & multiplier)
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

  bool XMLReader::ParseBool(std::string word)
  {
    bool ret = false;
    if( word == "False" || word == "false" || word == "0" || word == "no" || word == "N" || word == "No" || word == "NO" ) ret = false;
    else if( word == "True" || word == "true" || word == "1" || word == "yes" || word == "Y" || word == "Yes" || word == "YES" ) ret = true;
    else Report("Cannot understand boolean value %s",word.c_str());
    return ret;
  }

  void XMLReader::ParseCommaSeparated(std::string word, std::vector<std::string> & parsed, char symbol)
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

  void XMLReader::ParseReal(std::string word, std::vector<INMOST::Storage::real> & Vector, int & Repeat, int SetSize)
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

  void XMLReader::ParseVariable(std::string word, std::vector<INMOST::Storage::var> & Vector, int & Repeat, int SetSize)
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
  void XMLReader::ParseInteger(std::string word, std::vector<INMOST::Storage::integer> & Vector, int & Repeat, int SetSize)
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

  void XMLReader::ParseBulk(std::string word, std::string & Vector, int & Repeat, int SetSize)
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

  void XMLReader::ParseReference(std::string word, std::vector<std::pair<INMOST::ElementType,int> > & Vector, int & Repeat, int SetSize)
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

  void XMLReader::ParseRemoteReference(std::string word, std::vector< std::pair<std::string,std::pair<INMOST::ElementType,int> > > & Vector, int & Repeat, int SetSize)
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

  bool XMLReader::XMLTag::Failure() const {return finish == 0;}
  bool XMLReader::XMLTag::Process() const {return finish == 1 || finish == 2;}
  bool XMLReader::XMLTag::Finalize() const {return finish == 3;}
  const XMLReader::XMLAttrib & XMLReader::XMLTag::GetAttib(int n) const {return attributes[n];}
  XMLReader::XMLAttrib & XMLReader::XMLTag::GetAttib(int n) {return attributes[n];}
  int XMLReader::XMLTag::NumAttrib() const {return (int)attributes.size();}

  XMLReader::XMLTag XMLReader::OpenTag()
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

  bool XMLReader::CloseTag(XMLTag & tag)
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
}