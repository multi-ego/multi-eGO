import random

aforismi=[
    '“Imparerai a tue spese che nel lungo tragitto della vita incontrerai tante maschere e pochi volti.”',
    "“Mi si fissò invece il pensiero ch'io non ero per gli altri quel che finora, dentro di me, m'ero figurato d'essere.”",
    "“È molto più facileessere un eroe che un galantuomo. Eroi si può essere una volta tanto; galantuomini, si dev'esser sempre.”",
    "“Basta che lei si metta a gridare in faccia a tutti la verità. Nessuno ci crede, e tutti la prendono per pazza!”",
    "“La civiltà vuole che si auguri il buon giorno a uno che volentieri si manderebbe al diavolo; ed essere bene educati vuol dire appunto esser commedianti.”",
    "“Sapete che cosa significa amare l'umanità? Significa soltanto questo: essere contenti di noi stessi. Quando uno è contento di sé stesso, ama l'umanità.”",
    "“La facoltà d'illuderci che la realtà d'oggi sia la sola vera, se da un canto ci sostiene, dall'altro ci precipita in un vuoto senza fine, perché la realtà d'oggi é destinata a scoprire l'illusione domani. E la vita non conclude. Non può concludere. Se domani conclude, è finita.”",
    "“Ciò che conosciamo di noi è solamente una parte, e forse piccolissima, di ciò che siamo a nostra insaputa.”",
    "“Confidarsi con qualcuno, questo sì è veramente da pazzo!”",
    "“Se noi conosciamo che errare è dell'uomo non è crudeltà sovrumana la giustizia?”",
    "“Io non l'ho più questo bisogno, perché muojo ogni attimo io, e rinasco nuovo e senza ricordi: vivo e intero, non più in me, ma in ogni cosa fuori.”",
    "“Notiamo facilmente i difetti altrui e non ci accorgiamo dei nostri.”",
    "“Ogni cosa finché dura porta con sé la pena della sua forma, la pena d'esser così e di non poter essere più altrimenti.”",
]
    
def print_goodbye():
    print('-----------------------------===========+===+++++++++++++********#########################')
    print('-----------------------------=========++++**##########**********##########################')
    print('-----------------------------=-==+**#%%%%####*****#####%%%@%#########%%%%%################')
    print('----------------------=---====*#%@@@%%%%##**#******+****##%%@@@@@%%%%%%%##################')
    print('-------------------------==+#%@@@%%%###***+*+++=+*+*+++***##%@@@@@@@@@@%#%################')
    print('-----------------------==+%@@%%%###***+++========++=+++++**####%%@@@@@@@@%################')
    print('---------------------=+#%@@%###****+++++=============++=++++***###%@@@@@@@%###############')
    print('--------------------=*%%@@%###***+++====================+++++****###%@@@@@@%##############')
    print('------------------=*%@@@@%#***+++==========--============++++******##%%@@@@@@%############')
    print('-----------------=#@@@@@%#***++========-===--===--========+++++++****##%@@@@@@%#####%%%%%#')
    print('----------------+@@@@@%##**+++=====-==-=------------=======++++++++***##%%@@@@@###%%%%%%%%')
    print('---------------+@@@@@%#*+++=+=======--------==----------=========+**+*#%%%@@@@@%%%%%%%%%%%')
    print('--------------=%@@@@@##**++===-=---=--------------------=====++==+++***#%@@@@@@@%%%%@@@@@@')
    print('--------------+@@@@@%##*+++===----------------------------====+++++++*#%%%%@@@@@%%%%@@@@@@')
    print('-------------=#@@@@%%##+++===--------------------------------==+++++***#%%##@@@@%%%%%%@@@@')
    print('-------------=%@@@@%%#*+++=--------------------------------=====++++***#%%%%#@@@%%%%%%%%%%')
    print('-------------+@@@@@@###*++==-------------------------------====++++****#%%#%%%@@%#%%%%%%%%')
    print('------------=%@@@@@%###++===-=--------------------------======+++++****#%%%##%@@%######%%%')
    print('------------+@@@@@%%%#*++=====-==------------------=========++=*+++***####@@%%@@#########%')
    print('------------+@@@@@%%%**+======-=--------------------=======++++**********#%%%%@@########%%')
    print('------------+@@@@%%%##+++========------------------=-======+++++***+**###**###@@###*######')
    print('------------=%@@@%###*+=====------------------------======+++++++**++#*====+*#+=====+#####')
    print('------------=*@@@@#**+++==-===-------------=--------======++===+#+===+==----=+*+=---=*####')
    print('-------------=@@@@@@%*+====--=-------------------========+====+*=---==+*#*=--=+++=--=+%###')
    print('-------------+@@@@@@@@%#*=====---===----------------===++=====+==---#+=+*@#=--=+*+=--=#%##')
    print('-------------#@@@@@@@@@@@%#+======+======-========+*##%#+++=--===-=**#+=+======+++*--=%@##')
    print('------------=%@@@@@@@@@@@@@@%++===+=====++*#%%%@@@@%#*+=+%+=-------=#%+=====+=-===++-+#@@@')
    print('------------=@@@@@#@##@@@@@@@@#+=+*+==+*#%@@@%##*+=+**++#%+==-------=-==---==++===-=++##@@')
    print('------------=@%##%%**+###*#@@@@+=====+#%%%#%%@%%%%@@@@#*%++%#++---------=======*+==-+##%@@')
    print('------------+@##**#%%##**+++#@@+=--==++#%%%%%**@@@@%%%@@#**%@%+==---------===--=*+==+###@@')
    print('------------*@%###********++#@%==--===++=+*+==+++++=+++#%==+@#===---------===--=+++=+%%#@@')
    print('-----------=#@@##%###**###*#%%*=----=====*%%%#*++=====++@+==#*+---------==++===**#*+%@%%@@')
    print('-------===+#@@@%#*###**++##@@%+---====+--=+#%%%%##*#*++=#%+==---------==+***++*#%*##%@@@@@')
    print('---===*#%@@@@@@%*++======*%@@*=---=====-----==+*#*+=====+*=---------==+**+*+=+%%+**@@@@@@@')
    print('==+#%@@@@@@@@@@%*+======+*%@@*=---=====-----=======-===+===---------=++++++=+%#**#@@@@@@@@')
    print('#@@@@@@@@@@@@@@@#++=++*#%@%@@+=---====-==-=-====----=+*+==---------=====++=+#+*+#@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@#**##%%@#*@*==---=========+++====--++===--------=====++++=#+=+*@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@#*#@###==----====--==+=**+====+*===---------====+====*==*%@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@%+#@@@%+==---====--=+====*#*+++%+==---------==-===-=+#=+#@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@**+=+@@@@+=========+##+====+*%##@#+==------------==-=++=+#@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@%%#**%@@%*+*++#%@%#+=-=====++#@@#+==------------=-=*==-+#@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@##@@@@@@@@@%*==-------=-===%@#+===------------=+==-=+@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*===-=---======-+@@@#+=-------------=====+*%@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%%##+***%%*%%%%@@@@*=--------------==+=++@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%##**##*#+**+**+++*@@@@@@+--------------==+**%@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@%@@%#*#****#*+++*+=====+%@@@@@#=---------------=+**@%@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@*++#@%===-===============*@@@@@@+----------------==---*@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@*+=++===------========--=#%@@@@#=---------------------+@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@%=--=+==-------=---------*@@@@@%=----------------------*@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@*=-=------------------=+%@@@@@@=---------------------=*@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@#+=------------=====+%@@@@@@@@=------------------=*#@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@*=-----------==+*%@@@@@%#*++=----------------+#@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%=--------=*%%@@@@%#+==-----------------=+*%@@@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%+==--=+*%@@@@@@@=-------------------=*%@@@@@@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%##%@@@@@@@@@@*----------------=+#@@@@@@@@@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@+==-----==---=*#@@@@@@@@@@@@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@###*=---**+*%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%%%@@@@@@@@@@@@@@@@@#+*%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    print('')
    print('')
    print('')
    print('----------------------------------------------------')
    print('')
    print(':-)  LUIGI says goodbye  (-:')
    print('')
    print(aforismi[random.choice(range(len(aforismi)))])
    print()
    print('')
    print('----------------------------------------------------')
    print('')
    print('')
    print('')