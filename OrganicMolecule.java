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
                if (startChildren.contains(atom)) startChild = atom;
            }
            assignChildren(start, startChild);
        }
        
        private void assignChildren(IAtom a, IAtom parent)
        {
            
        }
        
        private IAtom atomWithLongestChain()
        {
            return null;
        }
        
        private boolean isNonAlkyl(IAtom a)
        {
            for (IBond bond : getConnectedBondsList(a))
            {
                for (IAtom atom : bond.atoms())
                {
                    if (atom.getAtomicNumber() != 6) return false;
                }
                if (bond.getOrder() != IBond.Order.SINGLE) return false;
            }
            return true;
        }
        
        public IAtom[] getChain()
        {
            IAtom[] chainArray = new IAtom[chain.size()];
            for (int i = 0; i < chainArray.length; i++)
                chainArray[i] = chain.get(i);
            return chainArray;
        }
    }
}