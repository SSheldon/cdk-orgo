import org.openscience.cdk.Molecule;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.exception.CDKException;
import java.util.*;

public class OrganicMolecule extends Molecule
{
    public enum Group { Alkyl, Alkenyl, Alkynyl, Imine, Amine, Sulfhydryl, Hydroxyl, Carbonyl, Formyl, Nitrile, Carboxyl }
    
    private IAtom[] chain;
    private String name;
    private Group principal;
    
    public OrganicMolecule(Molecule mol) throws CDKException
    {
        super(mol);
        boolean containsCarbon = false;
        for (IAtom atom : atoms())
        {
            if (atom.getSymbol().equals("C")) containsCarbon = true;
            break;
        }
        if (!containsCarbon) throw new CDKException("Molecule is inorganic.");
        chain = new CarbonChainFinder().getChain();
        principal = principalGroup();
        if (!verifyOrder())
        {
            for (int i = 0; i < chain.length / 2; i++)
            {
                IAtom temp = chain[i];
                chain[i] = chain[chain.length - 1 - i];
                chain[chain.length - 1 - i] = temp;
            }
        }
        name = new OrganicMoleculeNamer().getName();
    }
    
    public IAtom getChainCarbon(int position)
    { return chain[position]; }
    
    public int getChainLength()
    { return chain.length; }
    
    public String getIUPACName()
    { return name; }

    private boolean verifyOrder()
    {
        List<Integer> locs = principalGroupLocations();
        if (locs.size() == 0) return true;
        int first = locs.get(0);
        int last = locs.get(0);
        for (int i : locs)
        {
            first = Math.min(first, i);
            last = Math.max(last, i);
        }
        switch (((Integer)first).compareTo(chain.length - 1 - last))
        {
            case -1:
                return true;
            case 1:
                return false;
            case 0:
                int forwardSum = 0;
                int backwardSum = 0;
                for (int i : locs)
                {
                    forwardSum += i;
                    backwardSum += chain.length - 1 - i;
                }
                return forwardSum <= backwardSum;
            default:
                return first <= last;
        }
    }

    private boolean hasBond(int index, String symbol, int order)
    {
        for (IBond bond : getConnectedBondsList(chain[index]))
        {
            if (bond.getConnectedAtom(chain[index]).getSymbol().equals(symbol) && bond.getOrder().ordinal() + 1 == order) return true;
        }
        return false;
    }
    
    private Group highestPrecedenceGroup(int i)
    {
        if (!chain[i].getSymbol().equals("C")) throw new IllegalArgumentException();
        if (hasBond(i, "O", 2))
        {
            if (hasBond(i, "O", 1)) return Group.Carboxyl;
            else if (chain[i].getHydrogenCount() > 0) return Group.Formyl;
            else return Group.Carbonyl;
        }
        else if (hasBond(i, "N", 3)) return Group.Nitrile;
        else if (hasBond(i, "O", 1)) return Group.Hydroxyl;
        else if (hasBond(i, "S", 1)) return Group.Sulfhydryl;
        else if (hasBond(i, "N", 1)) return Group.Amine;
        else if (hasBond(i, "N", 2)) return Group.Imine;
        else if (hasBond(i, "C", 3)) return Group.Alkynyl;
        else if (hasBond(i, "C", 2)) return Group.Alkenyl;
        else return Group.Alkyl;
    }

    private Group principalGroup()
    {
        Group highest = Group.Alkyl;
        for (int i = 0; i < chain.length; i++)
        {
            Group g = highestPrecedenceGroup(i);
            if (g.compareTo(highest) > 0) highest = g;
        }
        return highest;
    }

    private List<Integer> principalGroupLocations()
    {            
        Group highest = principalGroup();
        if (highest == Group.Alkyl) return sideChainLocations();
        else
        {
            List<Integer> locs = new ArrayList<Integer>();
            for (int i = 0; i < chain.length; i++)
            {
                if (highestPrecedenceGroup(i) == highest) locs.add(i);
            }
            return locs;
        }
    }
    
    private List<Integer> sideChainLocations()
    {
        List<Integer> locs = new ArrayList<Integer>();
        for (int i = 0; i < chain.length; i++)
        {
            for (IBond bond : getConnectedBondsList(chain[i]))
            {
                if (bond.getOrder() == IBond.Order.SINGLE && !isInChain(bond.getConnectedAtom(chain[i]), i)) locs.add(i);
            }
        }
        return locs;
    }

    private boolean isInChain(IAtom b, int parentIndex)
    {
        if (!b.getSymbol().equals("C")) return false;
        if (parentIndex != 0 && b == chain[parentIndex - 1]) return true;
        if (parentIndex != chain.length - 1 && b == chain[parentIndex + 1]) return true;
        return false;
    }
    
    private class CarbonChainFinder
    {
        private List<IAtom> chain = new ArrayList<IAtom>();
        
        private void assignCarbonChain()
        {
            IAtom start = null;
            for (IAtom atom : atoms())
            {
                if (atom.getSymbol().equals("C") && isNonAlkyl(atom))
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
                if (child != parent && child.getSymbol().equals("C"))
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
                    if (!target.getSymbol().equals("C") || bond.getOrder() != IBond.Order.SINGLE) return true;
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
                if (child != parent && child.getSymbol().equals("C"))
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
                if (child != parent && child.getSymbol().equals("C"))
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
            for (IAtom a : atoms())
            {
                if (a.getSymbol().equals("C"))
                {
                    int length = longestChildChain(a, null);
                    if (length > maxLength)
                    {
                        maxLength = length;
                        atom = a;
                    }
                }
            }
            return atom;
        }
        
        private boolean isNonAlkyl(IAtom a)
        {
            if (!a.getSymbol().equals("C")) return true;
            for (IBond bond : getConnectedBondsList(a))
                if (!bond.getConnectedAtom(a).getSymbol().equals("C") || bond.getOrder() != IBond.Order.SINGLE) return true;
            return false;
        }
        
        private boolean verifyChain()
        {
            for (IAtom atom : atoms())
            {
                if (atom.getSymbol().equals("C") && isNonAlkyl(atom) && !chain.contains(atom)) return false;
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
        
        public IAtom[] getChain() throws CDKException
        {
            assignCarbonChain();
            if (!verifyChain()) throw new CDKException("Not all of the molecule's functional groups are included in the chain.");
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
                    if (a.getSymbol().equals("C") && a != parent && chain.contains(a))
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
    
    private class OrganicMoleculeNamer
    {
        private Map<String, List<Integer>> prefixes = new Hashtable<String, List<Integer>>();
        private Map<String, List<Integer>> suffixes = new Hashtable<String, List<Integer>>();

        private void enumerateChain()
        {
            for (int i = 0; i < chain.length; i++)
            {
                for (IBond bond : getConnectedBondsList(chain[i]))
                {
                    String symbol = bond.getConnectedAtom(chain[i]).getSymbol();
                    if (symbol.equals("O"))
                    {
                        if (bond.getOrder() == IBond.Order.DOUBLE)
                        {
                            if (hasBond(i, "O", 1)) addGroup(Group.Carboxyl, i);
                            else if (chain[i].getHydrogenCount() > 0) addGroup(Group.Formyl, i);
                            else addGroup(Group.Carbonyl, i);
                        }
                        else if (bond.getOrder() == IBond.Order.SINGLE && !hasBond(i, "O", 2)) addGroup(Group.Hydroxyl, i);
                    }
                    else if (symbol.equals("N"))
                    {
                        switch (bond.getOrder())
                        {
                            case SINGLE:
                                addGroup(Group.Amine, i);
                                break;
                            case DOUBLE:
                                addGroup(Group.Imine, i);
                                break;
                            case TRIPLE:
                                addGroup(Group.Nitrile, i);
                                break;
                        }
                    }
                    else if (symbol.equals("C"))
                    {
                        if (bond.getOrder() != IBond.Order.SINGLE)
                        {
                            if (i != chain.length - 1 && bond.getConnectedAtom(chain[i]) == chain[i + 1])
                            {
                                switch (bond.getOrder())
                                {
                                    case DOUBLE:
                                        addSuffix("en", i);
                                        break;
                                    case TRIPLE:
                                        addSuffix("yn", i);
                                        break;
                                }
                            }
                        }
                        else if (!isInChain(bond.getConnectedAtom(chain[i]), i))
                            addPrefix(chainLengthPrefix(alkylLength(bond.getConnectedAtom(chain[i]), chain[i])) + "yl", i);
                    }
                    else if (symbol.equals("F")) addPrefix("fluoro", i);
                    else if (symbol.equals("Cl")) addPrefix("chloro", i);
                    else if (symbol.equals("Br")) addPrefix("bromo", i);
                    else if (symbol.equals("I")) addPrefix("iodo", i);
                    else if (symbol.equals("S"))
                    {
                        if (bond.getOrder() == IBond.Order.SINGLE) addGroup(Group.Sulfhydryl, i);
                    }
                }
            }
        }

        private void addGroup(Group group, int index)
        {
            if (group == Group.Alkenyl || group == Group.Alkynyl || group == principal)
                addSuffix(groupSuffix(group), index);
            else
                addPrefix(groupPrefix(group), index);
        }

        private void addPrefix(String prefix, int index)
        {
            addAffix(prefixes, prefix, index);
        }

        private void addSuffix(String suffix, int index)
        {
            addAffix(suffixes, suffix, index);
        }

        private void addAffix(Map<String, List<Integer>> affixes, String affix, int index)
        {
            if (affixes.containsKey(affix)) affixes.get(affix).add(index);
            else
            {
                List<Integer> locs = new ArrayList<Integer>();
                locs.add(index);
                affixes.put(affix, locs);
            }
        }

        public String getName()
        {
            enumerateChain();
            String suffix = getSuffix();
            int firstNumber;
            for (firstNumber = 0; firstNumber < suffix.length() && !Character.isDigit(suffix.charAt(firstNumber)) && suffix.charAt(firstNumber) != '-'; firstNumber++) ;
            if (firstNumber == suffix.length()) firstNumber = 0;
            int followingLetter;
            for (followingLetter = firstNumber; followingLetter < suffix.length() && !Character.isLetter(suffix.charAt(followingLetter)); followingLetter++) ;
            String stem = suffix.substring(firstNumber, followingLetter) + chainLengthPrefix(chain.length) +
                (firstLetterIsConsonant(suffix) ? "a" : "") + suffix.substring(0, firstNumber) + suffix.substring(followingLetter);
            if (stem.charAt(0) == '-') stem = stem.substring(1);
            return Suffix(getPrefix(), stem);
        }

        private String getPrefix()
        {
            String prefix = "";
            while (prefixes.size() != 0)
            {
                String group = lowest(prefixes.keySet());
                prefix = Suffix(prefix, Affix(group, prefixes.get(group)));
                prefixes.remove(group);
            }
            return prefix;
        }

        private int alkylLength(IAtom b, IAtom parent)
        {
            for (IAtom child : getConnectedAtomsList(b))
            {
                if (child != parent) return alkylLength(child, b) + 1;
            }
            return 1;
        }

        private String getSuffix()
        {
            String suffix = "";
            if (suffixes.containsKey("en")) suffix = Suffix(suffix, Affix("en", suffixes.get("en")));
            if (suffixes.containsKey("yn"))
            {
                String alkynylSuffix = Affix("yn", suffixes.get("yn"));
                if (firstLetterIsConsonant(alkynylSuffix)) suffix += "e";
                suffix = Suffix(suffix, alkynylSuffix);
            }
            if (suffix.length() == 0) suffix = "an";
            switch (principal)
            {
                case Alkyl:
                case Alkenyl:
                case Alkynyl:
                    return suffix + "e";
                default:
                    String principalGroupSuffix = groupSuffix(principal);
                    String molGroupSuffix;
                    if (principal == Group.Carboxyl || principal == Group.Formyl || principal == Group.Nitrile)
                        molGroupSuffix = groupCountPrefix(suffixes.get(principalGroupSuffix).size()) + principalGroupSuffix;
                    else molGroupSuffix = Affix(principalGroupSuffix, suffixes.get(principalGroupSuffix));
                    if (firstLetterIsConsonant(molGroupSuffix)) suffix += "e";
                    return Suffix(suffix, molGroupSuffix);
            }
        }

        private boolean firstLetterIsConsonant(String s)
        {
            for (char c : s.toCharArray())
            {
                if (Character.isLetter(c)) return !(c == 'a' || c == 'e' || c == 'i' || c == 'o' || c == 'u' || c == 'y');
            }
            return false;
        }

        private String Suffix(String stem, String suffix)
        {
            if (stem.length() != 0 && suffix.length() != 0 && Character.isDigit(suffix.charAt(0))) return stem + "-" + suffix;
            else return stem + suffix;
        }

        private String prefix(String stem, String prefix)
        {
            if (stem.length() != 0 && Character.isDigit(stem.charAt(0))) return prefix + "-" + stem;
            else return prefix + stem;
        }

        private String Affix(String affix, List<Integer> locations)
        {
            Collections.sort(locations);
            String s = "";
            s += locations.get(0) + 1;
            for (int i = 1; i < locations.size(); i++)
                s += "," + (locations.get(i) + 1);
            s += "-" + groupCountPrefix(locations.size()) + affix;
            return s;
        }

        private String lowest(Set<String> affixes)
        {
            String lowest = null;
            for (String affix : affixes)
            {
                if (lowest == null || affix.compareTo(lowest) < 0) lowest = affix;
            }
            return lowest;
        }

        private String chainLengthPrefix(int digit, int place)
        {
            switch (place)
            {
                case 0:
                    switch (digit)
                    {
                        case 0: return "";
                        case 1: return "hen";
                        case 2: return "do";
                        case 3: return "tri";
                        case 4: return "tetr";
                        case 5: return "pent";
                        case 6: return "hex";
                        case 7: return "hept";
                        case 8: return "oct";
                        case 9: return "non";
                    }
                    throw new IllegalArgumentException();
                case 1:
                    switch (digit)
                    {
                        case 0: return "";
                        case 1: return "dec";
                        case 2: return "cos";
                        default: return chainLengthPrefix(digit, 0) + "acont";
                    }
                default: throw new IllegalArgumentException();
            }
        }

        private String chainLengthPrefix(int length)
        {
            switch (length)
            {
                case 1: return "meth";
                case 2: return "eth";
                case 3: return "prop";
                case 4: return "but";
                case 11:
                    return "un" + chainLengthPrefix(1, 1);
                case 20:
                    return "i" + chainLengthPrefix(2, 1);
                case 21:
                    return chainLengthPrefix(1, 0) + "i" + chainLengthPrefix(2, 1);
            }
            String s = "";
            for (int i = 0; length > 0; i++)
            {
                int digit = length % 10;
                if (digit != 0)
                {
                    s += chainLengthPrefix(digit, i);
                    if (i == 0 && length > 10 && digit > 3) s += "a";
                }
                length /= 10;
            }
            return s;
        }

        private String groupCountPrefix(int count)
        {
            switch (count)
            {
                case 1: return "";
                case 2: return "di";
                case 3: return "tri";
                case 4: return "tetra";
                case 5: return "penta";
            }
            throw new IllegalArgumentException();
        }

        private String groupPrefix(Group group)
        {
            switch (group)
            {
                case Imine: return "imino";
                case Amine: return "amino";
                case Sulfhydryl: return "sulfanyl";
                case Hydroxyl: return "hydroxy";
                case Carbonyl: return "oxo";
                case Formyl: return "formyl";
                case Nitrile: return "cyano";
                case Carboxyl: return "carboxy";
            }
            throw new IllegalArgumentException();
        }

        private String groupSuffix(Group group)
        {
            switch (group)
            {
                case Imine: return "imine";
                case Amine: return "amine";
                case Sulfhydryl: return "thiol";
                case Hydroxyl: return "ol";
                case Carbonyl: return "one";
                case Formyl: return "al";
                case Nitrile: return "nitrile";
                case Carboxyl: return "oic acid";
            }
            throw new IllegalArgumentException();
        }
    }
}