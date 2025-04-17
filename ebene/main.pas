unit main;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, math;

type
  TForm1 = class(TForm)
    GroupBox1: TGroupBox;
    Edit1: TEdit;
    Label1: TLabel;
    Label2: TLabel;
    Edit2: TEdit;
    Edit3: TEdit;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Edit4: TEdit;
    Button1: TButton;
    ListBox1: TListBox;
    GroupBox2: TGroupBox;
    Button2: TButton;
    Button3: TButton;
    OpenDialog1: TOpenDialog;
    SaveDialog1: TSaveDialog;
    Label8: TLabel;
    Button4: TButton;
    Button5: TButton;
    Button6: TButton;
    GroupBox3: TGroupBox;
    Memo1: TMemo;
    Button7: TButton;
    CheckBox1: TCheckBox;
    CheckBox2: TCheckBox;
    Edit5: TEdit;
    Button8: TButton;
    SaveDialog2: TSaveDialog;
    Button9: TButton;
    Button10: TButton;
    Button11: TButton;
    procedure Button1Click(Sender: TObject);
    procedure Button4Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure Button3Click(Sender: TObject);
    procedure Button7Click(Sender: TObject);
    procedure Button6Click(Sender: TObject);
    procedure Button5Click(Sender: TObject);
    procedure trafo(n:integer);
    procedure quadrate;
    procedure quadrate2;
    procedure Button8Click(Sender: TObject);
    procedure Button9Click(Sender: TObject);
    procedure Button10Click(Sender: TObject);
    procedure Button11Click(Sender: TObject);

  private
    { Private declarations }
  public
    { Public declarations }
  end;

  type Tgalaxien=record
  l,b,r:extended;
  name:string;
  end;

  type Tposi=record
  x,y,z:extended;
  end;

  type T2Darray=record
  data:array of extended;
  end;


var
  Form1: TForm1;
  galaxien:array of Tgalaxien;
  posi:array of Tposi;
  editon,sigmaclip,nurabstand:boolean;
  select:integer;
  nx,ny,nz,h,abw,abwmw,gb,gl:extended;

implementation

{$R *.dfm}

procedure TForm1.Button1Click(Sender: TObject);
begin
//liste erweitern   oder �ndern
if editon=false
then
begin
Listbox1.Items.add(edit4.Text);
setlength(galaxien,listbox1.Items.Count);
//Array f�llen
galaxien[listbox1.Items.Count-1].l:=strtofloat(edit1.text);
galaxien[listbox1.Items.Count-1].b:=strtofloat(edit2.text);
galaxien[listbox1.Items.Count-1].r:=strtofloat(edit3.text);
galaxien[listbox1.Items.Count-1].name:=edit4.text;
//Felder l�schen
end
else
begin
//Array f�llen
galaxien[select].l:=strtofloat(edit1.text);
galaxien[select].b:=strtofloat(edit2.text);
galaxien[select].r:=strtofloat(edit3.text);
galaxien[select].name:=edit4.text;


end;

editon:=false;
edit1.text:='';
edit2.text:='';
edit3.text:='';
edit4.text:='';
end;

procedure TForm1.Button4Click(Sender: TObject);
begin             //neue Liste
 if MessageDlg('Do you really want to create a new list? This command resets the whole window.',
    mtConfirmation, [mbYes, mbNo], 0) = mrYes then
    begin
    //alles zur�cksetzen
    memo1.Clear;
    listbox1.Clear;
    setlength(galaxien,0);
    edit1.text:='';
edit2.text:='';
edit3.text:='';
edit4.text:='';
edit5.Text:='3';
label8.caption:='name: -----';
checkbox1.Checked:=false;
checkbox2.Checked:=false;
nurabstand:=false;

editon:=false;
    end;

end;

procedure TForm1.Button2Click(Sender: TObject);
var a,i:integer;
begin
//markierten Eintrag aus List l�schen

editon:=false;
//suche Markierung
a:=-1;
for i:=0 to listbox1.Items.Count-1 do
begin
if listbox1.Selected[i]=true
then
a:=i;
end;

if a<>-1
then
begin
//aus liste l�schen
listbox1.Items.Delete(a);


//array verkleinern und daten verschieben
for i:=a+1 to Listbox1.Items.count do
galaxien[i-1]:=galaxien[i];

setlength(galaxien,listbox1.items.Count);

end;

end;



procedure TForm1.FormCreate(Sender: TObject);
begin                    //init
editon:=false;
nurabstand:=false;
end;

procedure TForm1.Button3Click(Sender: TObject);
var a,i:integer;
begin                     //Eintrag bearbeiten
editon:=false;

//suche Markierung
a:=-1;
for i:=0 to listbox1.Items.Count-1 do
begin
if listbox1.Selected[i]=true
then
a:=i;
end;

if a<>-1
then
begin           //hole daten aus array
editon:=true;
select:=a;
edit1.Text:=floattostr(galaxien[select].l);
edit2.Text:=floattostr(galaxien[select].b);
edit3.Text:=floattostr(galaxien[select].r);
edit4.Text:=galaxien[select].name;
end;

end;




procedure Tform1.trafo(n:integer);
begin                                  //Trafo von galaktischen Koord in kartesische

posi[n].x:=galaxien[n].r*sin((90-galaxien[n].b)*PI/180)*cos(galaxien[n].l*PI/180);

posi[n].y:=galaxien[n].r*sin((90-galaxien[n].b)*PI/180)*sin(galaxien[n].l*PI/180);

posi[n].z:=galaxien[n].r*cos((90-galaxien[n].b)*PI/180);

end;


procedure TForm1.quadrate;
var i,ii,iii:integer;
    MatrixB,LoesMatrix:array[0..2,0..2] of extended;
    MatrixA:array of T2Darray;
    VektorH:array of extended;
    VektorG,loesung:array[0..2] of extended;
    det,loesdet,norm,fehler,ddd:extended;
begin
Memo1.Clear;

setlength(posi,listbox1.Items.count);
for i:=0 to listbox1.Items.Count-1 do
trafo(i);




 //kleinste Quadrate


setlength(MatrixA,listbox1.items.count);
for i:=0 to listbox1.Items.Count-1 do
setlength(MatrixA[i].data,3);

setlength(VektorH,listbox1.items.count);

for i:=0 to listbox1.Items.count-1 do
begin
VektorH[i]:=posi[i].x;
MatrixA[i].data[0]:=1;
MatrixA[i].data[1]:=posi[i].y;
MatrixA[i].data[2]:=posi[i].z;
end;

for i:=0 to 2 do
for ii:=0 to 2 do
begin
MatrixB[i,ii]:=0;
for iii:=0 to listbox1.items.Count-1 do
begin
MatrixB[i,ii]:=MatrixB[i,ii]+MatrixA[iii].data[i]*MatrixA[iii].data[ii];
end;
end;

for i:=0 to 2 do
begin
VektorG[i]:=0;
for ii:=0 to listbox1.Items.Count-1 do
begin
VektorG[i]:=VektorG[i]+MatrixA[ii].data[i]*VektorH[ii];
end;
end;


det:=MatrixB[0,0]*MatrixB[1,1]*MatrixB[2,2];
det:=det+MatrixB[0,1]*MatrixB[1,2]*MatrixB[2,0];
det:=det+MatrixB[0,2]*MatrixB[1,0]*MatrixB[2,1];
det:=det-MatrixB[0,2]*MatrixB[1,1]*MatrixB[2,0];
det:=det-MatrixB[0,1]*MatrixB[1,0]*MatrixB[2,2];
det:=det-MatrixB[0,0]*MatrixB[1,2]*MatrixB[2,1];

for i:=0 to 2 do
begin

LoesMatrix:=MatrixB;

for ii:=0 to 2 do
LoesMatrix[i,ii]:=VektorG[ii];

loesdet:=LoesMatrix[0,0]*LoesMatrix[1,1]*LoesMatrix[2,2];
loesdet:=loesdet+LoesMatrix[0,1]*LoesMatrix[1,2]*LoesMatrix[2,0];
loesdet:=loesdet+LoesMatrix[0,2]*LoesMatrix[1,0]*LoesMatrix[2,1];
loesdet:=loesdet-LoesMatrix[0,2]*LoesMatrix[1,1]*LoesMatrix[2,0];
loesdet:=loesdet-LoesMatrix[0,1]*LoesMatrix[1,0]*LoesMatrix[2,2];
loesdet:=loesdet-LoesMatrix[0,0]*LoesMatrix[1,2]*LoesMatrix[2,1];

loesung[i]:=loesdet/det;


end;



norm:=sqrt(1+loesung[1]*loesung[1]+loesung[2]*loesung[2]);
nx:=1/norm;
ny:=-loesung[1]/norm;
nz:=-loesung[2]/norm;
h:=loesung[0]/norm;






gl:=180/PI*arctan2(-ny,-nx);
if gl<0
then
gl:=360+gl;

gb:=90-180/PI*arctan2(Sqrt(nx*nx+ny*ny),-nz);









//Fehlerrechnung
abw:=0;
for i:=0 to listbox1.Items.Count-1 do
begin
ddd:=nx*posi[i].x+ny*posi[i].y+nz*posi[i].z-h;

if nurabstand=true
then
begin
Memo1.Lines.Add(floattostr(ddd));
end;

abw:=abw+sqr(ddd);
end;

nurabstand:=false;

abw:=sqrt(abw/(listbox1.items.count-1));
abwmw:=abw/sqrt((listbox1.items.count-1));


if sigmaclip=true
then

begin
sigmaclip:=false;

iii:=listbox1.items.count-1;

for i:=iii downto 0 do
begin
fehler:=abs(nx*posi[i].x+ny*posi[i].y+nz*posi[i].z-h);

if fehler>(strtofloat(edit5.text)*abw)
then
begin

listbox1.Items.Delete(i);

//array verkleinern und daten verschieben
for ii:=i+1 to listbox1.Items.Count do
galaxien[ii-1]:=galaxien[ii];

setlength(galaxien,listbox1.items.Count);

end;

end;

quadrate;

end;




end;

















                    //HIER WEITERPROGRAMMIERN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
procedure TForm1.quadrate2;            //2 auf der Ebene
var i,ii,iii:integer;
    b,g,loesung:extended;
    VektorV,VektorU:array of extended;
    norm,fehler:extended;
    gal1,gal2:array[0..2] of extended;
    l:array[0..3] of extended;
begin
Memo1.Clear;

setlength(posi,listbox1.Items.count);
for i:=0 to listbox1.Items.Count-1 do
trafo(i);

 //kleinste Quadrate


Setlength(VektorV,listbox1.items.count-2);
Setlength(VektorU,listbox1.items.count-2);

gal1[0]:=posi[0].x;
gal1[1]:=posi[0].y;
gal1[2]:=posi[0].z;

gal2[0]:=posi[1].x;
gal2[1]:=posi[1].y;
gal2[2]:=posi[1].z;



for i:=2 to listbox1.Items.count-1 do
begin

VektorU[i-2]:=-gal1[1]*posi[i].x-gal1[0]*posi[i].y+gal1[1]*gal2[0];
VektorU[i-2]:=VektorU[i-2]+gal2[1]*posi[i].x+gal1[0]*posi[i].y-gal1[0]*gal2[1];

VektorV[i-2]:=-gal1[2]*posi[i].x-gal2[0]*posi[i].z+gal1[2]*gal2[0];
VektorV[i-2]:=VektorV[i-2]+gal2[2]*posi[i].x+gal1[0]*posi[i].z-gal1[0]*gal2[2];

end;

b:=0;
g:=0;
for i:=0 to listbox1.items.count-3 do
begin
b:=b+VektorV[i]*VektorV[i];
g:=g+VektorU[i]*VektorV[i];
end;

loesung:=g/b;

l[3]:=-loesung;
l[2]:=1;
l[1]:=-(l[2]*(gal1[1]-gal2[1])+l[3]*(gal1[2]-gal2[2]))/(gal1[0]-gal2[0]);

l[0]:=l[1]*gal1[0]+l[2]*gal1[1]+l[3]*gal1[2];

norm:=sqrt(l[1]*l[1]+l[2]*l[2]+l[3]*l[3]);

nx:=l[1]/norm;
ny:=l[2]/norm;
nz:=l[3]/norm;
h:=l[0]/norm;


gl:=180/PI*arctan2(ny,nx);

gb:=90-180/PI*arctan2(Sqrt(nx*nx+ny*ny),nz);


//Fehlerrechnung
abw:=0;
for i:=2 to listbox1.Items.Count-1 do
abw:=abw+sqr(nx*posi[i].x+ny*posi[i].y+nz*posi[i].z-h);

abw:=sqrt(abw/(listbox1.items.count-3));
abwmw:=abw/sqrt((listbox1.items.count-3));



if sigmaclip=true
then

begin
sigmaclip:=false;

iii:=listbox1.items.count-1;

for i:=iii downto 0 do
begin
fehler:=abs(nx*posi[i].x+ny*posi[i].y+nz*posi[i].z-h);

if fehler>strtofloat(edit5.text)*abw
then
begin

listbox1.Items.Delete(i);

//array verkleinern und daten verschieben
for ii:=i+1 to listbox1.Items.Count do
galaxien[ii-1]:=galaxien[ii];

setlength(galaxien,listbox1.items.Count);

end;

end;

quadrate2;

end;




end;













procedure TForm1.Button7Click(Sender: TObject);

begin                       //Berechnen

Memo1.Clear;

sigmaclip:=checkbox1.checked;

if checkbox2.checked=false
then
begin

quadrate;
Memo1.Lines.Add('number of datapoints: '+inttostr(listbox1.items.count));

end
else        //erste 2 Galaxien auf Ebene
begin

quadrate2;
Memo1.Lines.Add('number of datapoints: '+inttostr(listbox1.items.count-2));
end;




Memo1.Lines.Add('-------------------------------------------');
Memo1.Lines.Add('plane: nx*x+ny*y+nz*z=h');
Memo1.Lines.Add('nx = '+floattostr(nx));
Memo1.Lines.Add('ny = '+floattostr(ny));
Memo1.Lines.Add('nz = '+floattostr(nz));
Memo1.Lines.Add('h = '+floattostr(h));
Memo1.Lines.Add('.....');
Memo1.Lines.Add('galactic latitude = '+floattostr(gb));
Memo1.Lines.Add('galactic longitude = '+floattostr(gl));
Memo1.Lines.Add('.....');
Memo1.Lines.Add('sigma:'+floattostr(abw));
Memo1.Lines.Add('mean error:'+floattostr(abwmw));

end;

procedure TForm1.Button6Click(Sender: TObject);
var TF:TextFile;
    i:integer;
begin           //Liste speichern

if savedialog1.Execute
then
begin

assignfile(TF,savedialog1.FileName);
label8.caption:='name: '+savedialog1.filename;
rewrite(TF);
for i:=0 to listbox1.Items.Count-1 do
begin
writeln(TF,galaxien[i].name);
writeln(TF,floattostr(galaxien[i].l));
writeln(TF,floattostr(galaxien[i].b));
writeln(TF,floattostr(galaxien[i].r));
writeln(TF,'----------------------');
end;

closefile(TF);
end;

end;

procedure TForm1.Button5Click(Sender: TObject);
var TF:TextFile;
    i:integer;
    s:string;
begin           //Liste laden

if opendialog1.Execute
then
begin
                      //l�schen
    memo1.Clear;
    listbox1.Clear;
    setlength(galaxien,0);
    edit1.text:='';
edit2.text:='';
edit3.text:='';
edit4.text:='';
edit5.Text:='3';
label8.caption:='name: '+opendialog1.filename;
checkbox1.Checked:=false;
checkbox2.Checked:=false;
nurabstand:=false;

editon:=false;


                          //datei lesen
assignfile(TF,opendialog1.FileName);
reset(TF);
i:=0;
While EOF(TF)=false DO
begin
setlength(galaxien,i+1);
readln(TF,s);
galaxien[i].name:=s;
listbox1.Items.Add(s);

readln(TF,s);
galaxien[i].l:=strtofloat(s);
readln(TF,s);
galaxien[i].b:=strtofloat(s);
readln(TF,s);
galaxien[i].r:=strtofloat(s);
readln(TF,s);
inc(i);
end;


closefile(TF);


end;

end;



{abw:=0;
for i:=0 to listbox1.items.count-1 do
abw:=abw+sqr(-posi[i].x+loesung[0]+loesung[1]*posi[i].y+loesung[2]*posi[i].z);

abw:=sqrt(abw/(listbox1.items.count-1));

Q[0]:=MatrixB[1,1]*MatrixB[2,2]-MatrixB[1,2]*MatrixB[2,1];
Q[1]:=MatrixB[0,0]*MatrixB[2,2]-MatrixB[0,2]*MatrixB[2,0];
Q[2]:=MatrixB[1,1]*MatrixB[0,0]-MatrixB[1,0]*MatrixB[0,1];

for i:=0 to 2 do
begin
Q[i]:=abs(Q[i]/det);
Q[i]:=Sqrt(Q[i])*abw;
end;

fh:=Q[0]/norm;
fy:=Q[1]/norm;
fz:=Q[2]/norm;
fx:=Sqrt(Sqr(fh/h)+Sqr(fy/ny)+Sqr(fz/nz))*nx;   }



procedure TForm1.Button8Click(Sender: TObject);
begin                    //Ergebnisse speichern

if savedialog2.Execute
then
Memo1.Lines.SaveToFile(savedialog2.filename);

end;

procedure TForm1.Button9Click(Sender: TObject);
begin //Schlie�en
close;
end;

procedure TForm1.Button10Click(Sender: TObject);
var i:integer;
begin
Memo1.Clear;
Memo1.Lines.Add('x           y            z');
setlength(posi,listbox1.Items.count);
for i:=0 to listbox1.Items.Count-1 do
begin
trafo(i);
Memo1.Lines.Add(floattostr(posi[i].x)+' '+floattostr(posi[i].y)+' '+floattostr(posi[i].z));
end;

end;

procedure TForm1.Button11Click(Sender: TObject);
begin
Memo1.Clear;
nurabstand:=true;
quadrate;
end;

end.
