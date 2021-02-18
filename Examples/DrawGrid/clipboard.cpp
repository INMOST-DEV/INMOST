#include "clipboard.h"

#if defined(_WIN32)
#include <windows.h>
bool  setTextToPasteboard(std::string as)
{
    
    size_t reqLength = ::MultiByteToWideChar( CP_UTF8, 0, as.c_str(), (int)as.length(), 0, 0 );
    
    std::wstring text( reqLength, L'\0' );
    
    ::MultiByteToWideChar( CP_UTF8, 0, as.c_str(), (int)as.length(), &text[0], (int)text.length() );
    
    bool ok = false;
    if (OpenClipboard(NULL))
	{
        EmptyClipboard();
        HGLOBAL hClipboardData;
        size_t bytes = text.length()+1 * sizeof(wchar_t);
        
        hClipboardData = GlobalAlloc(GMEM_DDESHARE, bytes*2);
        wchar_t * pchData = (wchar_t*)GlobalLock(hClipboardData);
        
        wcscpy(pchData, text.c_str());
        
        GlobalUnlock(hClipboardData);
        SetClipboardData(CF_UNICODETEXT,hClipboardData);
        CloseClipboard();
        ok = true;
    }
    return ok;
}
std::string getTextFromPasteboard()
{
    std::string clipBoardText="";
    if (OpenClipboard(NULL))
	{
        HANDLE hClipboardData = GetClipboardData(CF_UNICODETEXT);
        if(IsClipboardFormatAvailable(CF_UNICODETEXT))
		{
            wchar_t * pszText =NULL;
            pszText = (wchar_t *)GlobalLock(hClipboardData);
            if (pszText == NULL)
			{
                
            }
			else
			{
                std::wstring  pchData = pszText;
                char * mbstr2 = new char[pchData.size()*4];
                size_t bytes = pchData.length()+1 * sizeof(wchar_t);
                WideCharToMultiByte(CP_UTF8,0,pchData.c_str(),bytes*2,mbstr2,bytes*2,NULL,NULL);
                clipBoardText.append(mbstr2);
				delete [] mbstr2;
            }
            GlobalUnlock(hClipboardData);
            CloseClipboard();
        }
    }
    return clipBoardText;
}
#elif defined(__APPLE__)

#include <Carbon/Carbon.h>
//AUTOFRAMEWORK(Carbon)


bool setTextToPasteboard(std::string str) 
{
    const char * byteArrayIndex = str.c_str();
    OSStatus err = noErr;
    static PasteboardRef	pasteboard = NULL;
    PasteboardCreate( kPasteboardClipboard, &pasteboard );
    
    err = PasteboardClear( pasteboard );
    //require_noerr( err, PasteboardClear_FAILED );
    
    CFDataRef data;
    
    data =   CFDataCreate(kCFAllocatorDefault, (UInt8*)byteArrayIndex, strlen(byteArrayIndex)+1);
    
    err = PasteboardPutItemFlavor( pasteboard, (PasteboardItemID)1, kUTTypeUTF8PlainText, data, 0);
    //require_noerr( err, PasteboardPutItemFlavor_FAILED );
    
PasteboardPutItemFlavor_FAILED:
PasteboardClear_FAILED:
    return err == noErr;
}
std::string getTextFromPasteboard() 
{
    
    std::string clipBoard = "";
    OSStatus err = noErr;
    ItemCount  itemCount;
    PasteboardSyncFlags  syncFlags;
    static PasteboardRef inPasteboard = NULL;
    PasteboardCreate( kPasteboardClipboard, &inPasteboard );
    char* data = NULL;
    
    syncFlags = PasteboardSynchronize( inPasteboard );
    err = badPasteboardSyncErr;
    
    err = PasteboardGetItemCount( inPasteboard, &itemCount );
    //require_noerr( err, CantGetPasteboardItemCount );
    
    for( int itemIndex = 1; itemIndex <= itemCount; itemIndex++ ) 
	{
        PasteboardItemID itemID;
        CFDataRef flavorData;
        err = PasteboardGetItemIdentifier( inPasteboard, itemIndex, &itemID );
        //require_noerr( err, CantGetPasteboardItemIdentifier );
        err = PasteboardCopyItemFlavorData( inPasteboard, itemID, CFSTR("public.utf8-plain-text"), &flavorData );
        if(err==noErr)data = (char*)CFDataGetBytePtr(flavorData);
        if( data!=NULL && err==noErr )
			clipBoard.append(data);
		else 
			return "Error Pasting";
CantGetPasteboardItemIdentifier:
        ;
    }
CantGetPasteboardItemCount:
    return clipBoard;
    
}
#else
bool setTextToPasteboard(std::string str) {(void)str; return false;}
std::string getTextFromPasteboard() {return "";}
#endif
