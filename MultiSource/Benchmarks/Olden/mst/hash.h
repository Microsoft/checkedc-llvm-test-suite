/* For copyright information, see olden_v1.0/COPYRIGHT */

#include <stdchecked.h>
#include "stdio.h"
#include <stdio_checked.h>

struct hash_entry {
  unsigned int key;
  void *entry;
  ptr<struct hash_entry> next;
};

typedef ptr<struct hash_entry> HashEntry;

struct hash {
  HashEntry *array : itype(array_ptr<HashEntry>);
  ptr<int(unsigned int)> mapfunc;
  int size;
};

typedef ptr<struct hash> Hash;

Hash MakeHash(int size, ptr<int(unsigned int)> map);
void *HashLookup(unsigned int key, Hash hash);
void HashInsert(void *entry, unsigned int key, Hash hash);
void HashDelete(unsigned int key, Hash hash);
