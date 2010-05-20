import org.openscience.cdk.Molecule;
import org.openscience.cdk.interfaces.*;
import java.util.*;

public class OrganicMolecule extends Molecule
{
    private IAtom[] chain;
    private String name;
    
    public OrganicMolecule(Molecule mol)
    {
        super(mol);
        chain = new CarbonChainFinder().getChain();
    }
    
    public IAtom getChainCarbon(int position)
    { return chain[position]; }
    
    public int getChainLength()
    { return chain.length; }
    
    public String getIUPACName()
    { return name; }
    
    private class CarbonChainFinder
    {
        private List<IAtom> chain = new ArrayList<IAtom>();
        
        private void assignCarbonChain()
        {
            IAtom start = null;
            for (IAtom atom : atoms())
            {
                if (atom.getAtomicNumber() == 6 && isNonAlkyl(atom))
                    start = atom;
            }
            if (start == null) start = atomWithLongestChain();
            chain.add(start);
            assignChildren(start, null);
            IAtom startChild = null;
            List<IAtom> startChildren = getConnectedAtomsList(start);
            for (IAtom atom : chain)
            {
                if (startChildren.contains(atom))
                {
                    startChild = atom;
                    break;
                }
            }
            startChildren = null;
            assignChildren(start, startChild);
        }
        
        private void assignChildren(IAtom a, IAtom parent)
        {
            for (IAtom child : getConnectedAtomsList(a))
            {
                if (child != parent && child.getAtomicNumber() == 6)
                {
                    if (hasNonAlkylChild(child, a))
                    {
                        chain.add(child);
                        assignChildren(child, a);
                        return;
                    }
                }
            }
            IAtom child = childWithLongestChain(a, parent);
            if (child != null)
            {
                chain.add(child);
                assignChildren(child, a);
            }
        }
        
        private boolean hasNonAlkylChild(IAtom a, IAtom parent)
        {
            for (IBond bond : getConnectedBondsList(a))
            {
                IAtom target = bond.getConnectedAtom(a);
                if (target != parent)
                {
                    if (target.getAtomicNumber() != 6 || bond.getOrder() != IBond.Order.SINGLE) return true;
                    if (hasNonAlkylChild(target, a)) return true;
                }
            }
            return false;
        }

        private IAtom childWithLongestChain(IAtom a, IAtom parent)
        {
            int max = 0;
            IAtom maxChild = null;
            for (IAtom child : getConnectedAtomsList(a))
            {
                if (child != parent && child.getAtomicNumber() == 6)
                {
                    int childChainLength = longestChildChain(child, a) + 1;
                    if (childChainLength > max) maxChild = child;
                    max = Math.max(max, childChainLength);
                }
            }
            return maxChild;
        }

        private int longestChildChain(IAtom a,IAtom parent)
        {
            int max = 0;
            for (IAtom child : getConnectedAtomsList(a))
            {
                if (child != parent && child.getAtomicNumber() == 6)
                {
                    max = Math.max(max, longestChildChain(child, a) + 1);
                }
            }
            return max;
        }
        
        private IAtom atomWithLongestChain()
        {
            int maxLength = 0;
            IAtom atom = null;
            for (IAtom a : chain)
            {
                int length = longestChildChain(a, null);
                if (length > maxLength)
                {
                    maxLength = length;
                    atom = a;
                }
            }
            return atom;
        }
        
        private boolean isNonAlkyl(IAtom a)
        {
            if (a.getAtomicNumber() != 6) return false;
            for (IBond bond : getConnectedBondsList(a))
            {
                if (bond.getConnectedAtom(a).getAtomicNumber() != 6) return false;
                if (bond.getOrder() != IBond.Order.SINGLE) return false;
            }
            return true;
        }
        
        private boolean verifyChain()
        {
            for (IAtom atom : atoms())
            {
                if (atom.getAtomicNumber() == 6 && isNonAlkyl(atom) && !chain.contains(atom)) return false;
            }
            return true;
        }

        private boolean isEnd(IAtom a)
        {
            int bondedCInChain = 0;
            for (IAtom atom : getConnectedAtomsList(a))
            {
                if (chain.contains(atom)) bondedCInChain++;
            }
            return bondedCInChain < 2;
        }
        
        public IAtom[] getChain()
        {
            assignCarbonChain();
            IAtom[] chainArray = new IAtom[chain.size()];
            for (IAtom a : chain)
            {
                if (isEnd(a))
                {
                    chainArray[0] = a;
                    break;
                }
            }
            IAtom parent = null;            
            for (int i = 1; i < chainArray.length; i++)
            {
                for (IAtom a : getConnectedAtomsList(chainArray[i - 1]))
                {
                    if (a.getAtomicNumber() == 6 && a != parent && chain.contains(a))
                    {
                        parent = chainArray[i - 1];
                        chainArray[i] = a;
                        break;
                    }
                }
            }
            return chainArray;
        }
    }
}