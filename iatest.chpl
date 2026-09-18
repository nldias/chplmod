use IO;
const fou = openWriter("iatest.out");
for i  in 1..10 do {
  fou.writef("%i\n",i);
}
fou.close();