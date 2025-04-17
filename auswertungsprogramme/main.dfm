object Form1: TForm1
  Left = 255
  Top = 180
  Width = 857
  Height = 598
  Caption = 'Form1'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Memo1: TMemo
    Left = 8
    Top = 8
    Width = 833
    Height = 497
    ScrollBars = ssBoth
    TabOrder = 0
  end
  object Button1: TButton
    Left = 8
    Top = 512
    Width = 353
    Height = 49
    Cursor = crHandPoint
    Caption = 'lade'
    TabOrder = 1
    OnClick = Button1Click
  end
  object Button2: TButton
    Left = 368
    Top = 512
    Width = 353
    Height = 49
    Cursor = crHandPoint
    Caption = 'speichern'
    TabOrder = 2
    OnClick = Button2Click
  end
  object SpinEdit1: TSpinEdit
    Left = 728
    Top = 512
    Width = 41
    Height = 22
    MaxValue = 45
    MinValue = 0
    TabOrder = 3
    Value = 41
  end
  object SpinEdit2: TSpinEdit
    Left = 728
    Top = 536
    Width = 41
    Height = 22
    MaxValue = 5
    MinValue = 0
    TabOrder = 4
    Value = 5
  end
  object Edit1: TEdit
    Left = 776
    Top = 512
    Width = 65
    Height = 21
    TabOrder = 5
    Text = '3'
  end
  object Edit2: TEdit
    Left = 776
    Top = 536
    Width = 65
    Height = 21
    TabOrder = 6
    Text = '378,669'
  end
  object OpenDialog1: TOpenDialog
    Left = 48
    Top = 64
  end
  object SaveDialog1: TSaveDialog
    Left = 80
    Top = 64
  end
end
