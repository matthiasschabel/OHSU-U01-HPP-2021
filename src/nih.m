% NIH color map with "nih" colormap from OsiriX (http://www.osirix-viewer.com)
%   NIH(M) returns an M-by-3 matrix containing a "nih" colormap.
%   NIH, by itself, is the same length as the current colormap.
function cmap = nih(m)
    persistent rgb;

    if (nargin < 1)
       m = size(get(gcf,'colormap'),1);
    end

    if (isempty(rgb))
        rgb = ...
        [    0     0     0;...
             5     0    11;...
            11     0    22;...
            17     0    34;...
            22     0    45;...
            28     0    56;...
            34     0    68;...
            39     0    79;...
            45     0    90;...
            51     0   102;...
            56     0   113;...
            62     0   124;...
            68     0   136;...
            73     0   147;...
            79     0   158;...
            85     0   170;...
            79     0   164;...
            74     0   159;...
            69     0   154;...
            63     0   148;...
            58     0   143;...
            53     0   138;...
            47     0   132;...
            42     0   127;...
            37     0   122;...
            31     0   116;...
            26     0   111;...
            21     0   106;...
            15     0   100;...
            10     0    95;...
             5     0    90;...
             0     0    85;...
             0     0    90;...
             0     0    95;...
             0     0   100;...
             0     0   106;...
             0     0   111;...
             0     0   116;...
             0     0   122;...
             0     0   127;...
             0     0   132;...
             0     0   138;...
             0     0   143;...
             0     0   148;...
             0     0   154;...
             0     0   159;...
             0     0   164;...
             0     0   170;...
             0     0   175;...
             0     0   180;...
             0     0   185;...
             0     0   191;...
             0     0   196;...
             0     0   201;...
             0     0   207;...
             0     0   212;...
             0     0   217;...
             0     0   223;...
             0     0   228;...
             0     0   233;...
             0     0   239;...
             0     0   244;...
             0     0   249;...
             0     0   255;...
             0     5   255;...
             0    10   255;...
             0    15   255;...
             0    21   255;...
             0    26   255;...
             0    31   255;...
             0    37   255;...
             0    42   255;...
             0    47   255;...
             0    53   255;...
             0    58   255;...
             0    63   255;...
             0    69   255;...
             0    74   255;...
             0    79   255;...
             0    85   255;...
             0    90   249;...
             0    95   244;...
             0   100   239;...
             0   106   233;...
             0   111   228;...
             0   116   223;...
             0   122   217;...
             0   127   212;...
             0   132   207;...
             0   138   201;...
             0   143   196;...
             0   148   191;...
             0   154   185;...
             0   159   180;...
             0   164   175;...
             0   170   170;...
             0   175   170;...
             0   180   170;...
             0   185   170;...
             0   191   170;...
             0   196   170;...
             0   201   170;...
             0   207   170;...
             0   212   170;...
             0   217   170;...
             0   223   170;...
             0   228   170;...
             0   233   170;...
             0   239   170;...
             0   244   170;...
             0   249   170;...
             0   255   170;...
             0   255   159;...
             0   255   148;...
             0   255   138;...
             0   255   127;...
             0   255   116;...
             0   255   106;...
             0   255    95;...
             0   255    85;...
             0   255    74;...
             0   255    63;...
             0   255    53;...
             0   255    42;...
             0   255    31;...
             0   255    21;...
             0   255    10;...
             0   255     0;...
             5   255     5;...
            10   255    10;...
            15   255    15;...
            21   255    21;...
            26   255    26;...
            31   255    31;...
            37   255    37;...
            42   255    42;...
            47   255    47;...
            53   255    53;...
            58   255    58;...
            63   255    63;...
            69   255    69;...
            74   255    74;...
            79   255    79;...
            85   255    85;...
            95   255    79;...
           106   255    74;...
           116   255    69;...
           127   255    63;...
           138   255    58;...
           148   255    53;...
           159   255    47;...
           170   255    42;...
           180   255    37;...
           191   255    31;...
           201   255    26;...
           212   255    21;...
           223   255    15;...
           233   255    10;...
           244   255     5;...
           255   255     0;...
           255   249     0;...
           255   244     0;...
           255   239     0;...
           255   233     0;...
           255   228     0;...
           255   223     0;...
           255   217     0;...
           255   212     0;...
           255   207     0;...
           255   201     0;...
           255   196     0;...
           255   191     0;...
           255   185     0;...
           255   180     0;...
           255   175     0;...
           255   170     0;...
           255   164     0;...
           255   159     0;...
           255   154     0;...
           255   148     0;...
           255   143     0;...
           255   138     0;...
           255   132     0;...
           255   127     0;...
           255   122     0;...
           255   116     0;...
           255   111     0;...
           255   106     0;...
           255   100     0;...
           255    95     0;...
           255    90     0;...
           255    85     0;...
           255    81     0;...
           255    78     0;...
           255    75     0;...
           255    71     0;...
           255    68     0;...
           255    65     0;...
           255    62     0;...
           255    58     0;...
           255    55     0;...
           255    52     0;...
           255    49     0;...
           255    45     0;...
           255    42     0;...
           255    39     0;...
           255    35     0;...
           255    32     0;...
           255    29     0;...
           255    26     0;...
           255    22     0;...
           255    19     0;...
           255    16     0;...
           255    13     0;...
           255     9     0;...
           255     6     0;...
           255     3     0;...
           255     0     0;...
           252     0     0;...
           250     0     0;...
           248     0     0;...
           246     0     0;...
           243     0     0;...
           241     0     0;...
           239     0     0;...
           237     0     0;...
           234     0     0;...
           232     0     0;...
           230     0     0;...
           228     0     0;...
           225     0     0;...
           223     0     0;...
           221     0     0;...
           219     0     0;...
           216     0     0;...
           214     0     0;...
           212     0     0;...
           210     0     0;...
           208     0     0;...
           205     0     0;...
           203     0     0;...
           201     0     0;...
           199     0     0;...
           196     0     0;...
           194     0     0;...
           192     0     0;...
           190     0     0;...
           187     0     0;...
           185     0     0;...
           183     0     0;...
           181     0     0;...
           178     0     0;...
           176     0     0;...
           174     0     0;...
           255   255   255;...
           255   255   255 ];
    end;

    cmap = zeros([m 3]);

    % red 
    cmap(:,1) = interp1(0:1/255:1,rgb(:,1),0:1/(m-1):1)/255;

    % green
    cmap(:,2) = interp1(0:1/255:1,rgb(:,2),0:1/(m-1):1)/255;

    % blue
    cmap(:,3) = interp1(0:1/255:1,rgb(:,3),0:1/(m-1):1)/255;
end       