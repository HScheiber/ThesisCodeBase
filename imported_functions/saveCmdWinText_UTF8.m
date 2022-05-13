function saveCmdWinText_UTF8(filename)
    cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
    txt = char(cmdWinDoc.getText(0,cmdWinDoc.getLength));
    fid = fopen(filename,'W','n','utf-8');
    fwrite(fid,txt,'char');
    fclose(fid);
    %winopen(filename);  % in case you wish to verify...
end