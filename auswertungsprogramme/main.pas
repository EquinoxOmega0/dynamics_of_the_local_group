unit main;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, Spin, StdCtrls, math;

type
  TForm1 = class(TForm)
    Memo1: TMemo;
    Button1: TButton;
    Button2: TButton;
    OpenDialog1: TOpenDialog;
    SaveDialog1: TSaveDialog;
    SpinEdit1: TSpinEdit;
    SpinEdit2: TSpinEdit;
    Edit1: TEdit;
    Edit2: TEdit;
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form1: TForm1;
  sv,sp:array[1..5,1..45] of string;
  s:string;
  galnam:array[1..45] of string;
  hv,hp,scalev,scalep:extended;
  F:textfile;

implementation

{$R *.dfm}

procedure TForm1.Button1Click(Sender: TObject);
var i,ii:integer;
begin

assignfile(F,'galname.txt');
Reset(F);
For i:=1 to 45 do
begin
Readln(F,galnam[i]);
end;
CloseFile(F);

if OpenDialog1.Execute then
begin
assignfile(F,opendialog1.FileName);
Reset(F);
For ii:=1 to spinedit2.value do
begin
For i:=1 to spinedit1.value do
begin
Readln(F,sp[ii,i]);
Readln(F,sv[ii,i]);
end;
Readln(F,s);
end;


closefile(F);

scalep:=strtofloat(edit1.text);
scalev:=strtofloat(edit2.text);

For i:=1 to spinedit1.value do
begin
s:='';
s:='\hline '+galnam[i];
for ii:=1 to spinedit2.value do
begin
hp:=strtofloat(sp[ii,i]);
hv:=strtofloat(sv[ii,i]);
hp:=hp*scalep;
hv:=hv*scalev;
hp:=round(hp);
hv:=round(hv);
s:=s+' & '+floattostr(hp)+' & '+floattostr(hv);
end;

s:=s+'\\';

Memo1.Lines.Add(s);
end;

end;

end;

procedure TForm1.Button2Click(Sender: TObject);
begin
if SaveDialog1.Execute then
begin
Memo1.Lines.SaveToFile(SaveDialog1.filename);
end;
end;

end.
