vmm_acq_online() {
    outfile="$(date +%F_%H-%M-%S).pcapng"
    echo "📁 Fichier : $outfile"
    filter="dst host 10.0.0.3"

    # Lancer dumpcap en foreground, Python dans un autre terminal ou via screen/tmux
    dumpcap -i enp2s0 -f "$filter" -w "$outfile" &
    DUMPCAP_PID=$!

    # Ctrl+C dans le shell arrête Python, ensuite on kill dumpcap
    trap "echo '🛑 Arrêt demandé'; kill $DUMPCAP_PID; wait $DUMPCAP_PID; exit" SIGINT SIGTERM

    python3 ~/scripts/FEC_Datadescrambler_Default.py "$outfile"
}