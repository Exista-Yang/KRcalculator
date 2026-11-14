import json,os,sys
# validation (b): C2-4 binary vs C2-2 results, closure by closure
tags=[l.split(';')[0] for l in open('runs_c2/k4/markov_words.txt') if l.strip()]
ok=0; n=0
for t in tags:
    for F in ['31','32b']:
        a='runs_c2/k4m/markov/%s_f%s.json'%(t,F); b='runs_c2/k4/markov/%s_f%s.json'%(t,F)
        if not os.path.exists(a): continue
        A=json.load(open(a))['results']; B=json.load(open(b))['results']
        same=len(A)==len(B) and all(x['N']==y['N'] and x['reduced']==y['reduced'] and x['P']['str']==y['P']['str'] for x,y in zip(A,B))
        n+=1; ok+=same
        print(t,F,'IDENTICAL' if same else 'DIFFER', [x['P']['rank'] for x in A])
print('%d/%d identical'%(ok,n))
